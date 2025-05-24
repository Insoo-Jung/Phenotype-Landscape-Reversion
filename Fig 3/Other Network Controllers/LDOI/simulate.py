'''
Simulate organizes dynamics of networks into steady states, including attractors and phenotypes
and calculates their basin sizes by exhaustive simulation	

Note that phenotypes are considered to be the states of input and output nodes, which are specified in a settings.yaml file					
last update 7/17/23	

'''

import itertools, util, math, sys
import lap, param, features 
from net import Net
from copy import deepcopy
from timeit import default_timer as timer
CUPY, cp = util.import_cp_or_np(try_cupy=1) #should import numpy as cp if cupy not installed
# anything going into lap.py (i.e. x0) should be cp, but other matrices should be np
import numpy as np


#########################################################################################################


def init(param_file):
	params = param.load(param_file)
	G = Net(params)
	return params, G


#########################################################################################################

class Attractor:
	def __init__(self, params, G, attr_data):
		# attr_data mandatory keys: id, state, size, avg

		# add features list, but keep id,state,size,avg for compatibility

		self.id = attr_data['id']		# note id is just the string version of state
		self.state = attr_data['state'] # this is ONE set in the attractor (even if the attractor actually oscilaltes)

		self.size = attr_data['size']
		if 'avg' in attr_data.keys():
			self.avg = attr_data['avg']

		self.features = attr_data

		if util.istrue(params,'use_phenos'):
			self.map_to_pheno(params,G) 

	def map_to_pheno(self, params, G):
		outputs = G.output_indices()

		if 'inputs' in params.keys():
			inputs = G.input_indices()

		self.phenotype = ''
		if 'inputs' in params.keys() and params['use_inputs_in_pheno']:
			thresh=.5
			for i in range(len(inputs)):
				if self.avg[inputs[i]] > thresh: 
					self.phenotype +='1'
				else:
					self.phenotype +='0'
			if len(inputs) > 0:	
				self.phenotype +='|'

		self.pheno_arr = []
		for i in range(len(outputs)):

			assert((len(outputs)==len(params['output_thresholds'])) | (len(params['output_thresholds'])==1)) # otherwise number of outputs != number of output thresholds!
			#if float(self.avg[outputs[i]]) > params['output_thresholds'][i]: 
			if len(params['output_thresholds']) == 1: 	# such that if only 1 output threshold is given, it is always used
				thresh_index=0
			else: 										# else one threshold per output
				thresh_index=i 
			if float(self.avg[outputs[i]]) > params['output_thresholds'][thresh_index]:
				self.phenotype +='1'
				self.pheno_arr += [1]
			else:
				self.phenotype +='0'
				self.pheno_arr += [0]
				
	
	def __str__(self):		 
		return str(self.id) # note that id is only ONE state in the attractor

class Phenotype:
	# Phenotypes are defined by both the input and output node states
	def __init__(self, pheno_id, attractors, size, inputs, outputs):
		self.id = pheno_id
		self.attractors = attractors # []
		self.size = size
		self.inputs = inputs 
		self.outputs = outputs

	def __str__(self):
		return 'Inputs ' + str(self.inputs) + '\tOutputs ' + str(self.outputs) + '\tsize ' + str(self.size)

		
class SteadyStates:
	# collects sets of attractors and phenotypes

	def __init__(self, params,G):
		self.attractors = [] 
		self.phenotypes = []
		self.stats = {}
		self.params = params
		self.G = G
		self.attractor_features = {}

		if util.istrue(params,'track_lap_stats'):
			self.stats = {}

	def print_phenos(self):
		print('Inputs',self.params['inputs'],'\nOutputs',self.params['outputs'])
		for i in range(len(self.phenotypes)):
			print(self.phenotypes[i])

	def phenos_str(self):
		assert(0) #need to fix
		s=''
		for k in self.phenotypes:
			if s!='':
				s+=', '
			s+=k+':'+str(self.phenotypes[k].size)
		return s

	def build_attractors(self, feats):
		# note that attractor features are stored both in SteadyState.attractor_features[feature_id][attractor_number]
		# and in each Attractor.features[feature_id], this later one being used less (mainly for back compat)

		vals = feats['attractor_id'].val
		if CUPY:
			vals = vals.get()
		#cupy_to_numpy(self.params, feats['attractor_id'])
		unique, indx, inverse, counts = np.unique(vals,axis=0,return_index=True,return_inverse=True, return_counts=True)
		# usage:
		# 	unique, indx, inverse, counts = np.unique(x)
		# 	x[indx] == unique
		# 	unique[inverse] = x
		#	note this function is not yet implemented in cupy


		# TODO: CLEAN THIS
		#	Q is if there will only be adding and bool type, or could later develop weirder features
		#		leaning twds just + and bool for now
		#		...but even bool, not all will use 'any', maybe specify a combining function in [all, any, mean, sum...]

		# jp want to include this attractor_ids string too as a 'feature'
		num_attractors = len(unique)
		attractor_ids = [format_id_str(unique[i]) for i in range(num_attractors)]	

		# should probably mv this into features.py, not sure exactly which parts yet tho
		for feature in feats.values():
			if util.istrue(feature.args,'attractor'):
				if feature.id not in ['attractor_id','id','state','size']: # handled seperately below
					#print('\nfeature id',feature.id,feature.val.shape,self.params['num_samples'],'\n\n')
					assert(len(feature.val) == self.params['num_samples']) # assumes first index is for samples
					#if feature.val.dtype == bool: # this assumes all features are np/cp
					#	dt = bool
					#else:
					dt = float # since avg bool too
					shape = (num_attractors,) + feature.val[0].shape # assumes 0th axis of feature is for samples (which are now being mapped to attractors)
					feature.attractor_val = np.zeros(shape, dtype=dt)
					self.attractor_features[feature.id] = np.zeros(shape, dtype=dt)

		# default attractor features
		self.attractor_features['id'] = [None for i in range(num_attractors)] # string ids
		self.attractor_features['state'] = np.zeros((num_attractors,self.G.n),dtype=bool) # a boolean state in the attractor, string version is id
		self.attractor_features['size'] = np.zeros(num_attractors) # fraction of states that end in this attractor

		for i in range(num_attractors):  #expect this to be << num_samples
			attr_features = {}
			# putting things in: 1) features, 2) attractors, and 3) attractor_features...
			for feature in feats.values():
				if feature.id not in ['attractor_id','id','state','size']: # handled seperately below
					if util.istrue(feature.args,'attractor'):
						#if feature.val.dtype == bool: # this assumes all features are np/cp
						#	feature.attractor_val[i] = np.any(feature.val[inverse==i], axis=0)
						#else:
						if CUPY:
							feature.attractor_val[i] = np.sum(feature.val.get()[inverse==i], axis=0)
						else:
							feature.attractor_val[i] = np.sum(feature.val[inverse==i], axis=0)
						feature.attractor_val[i] /= counts[i]
						attr_features[feature.id] = feature.attractor_val[i]
						self.attractor_features[feature.id][i] = feature.attractor_val[i]

			# add default features
			self.attractor_features['id'][i] = attractor_ids[i]
			self.attractor_features['state'][i] = unique[i]
			self.attractor_features['size'][i] = counts[i]/self.params['num_samples']
		
			# since Attractor objs break the cp form, may be applications where don't want it...
			# this time this will all be cp except poss id
			attr_data = {'id':attractor_ids[i], 'state':unique[i], 'size':counts[i]/self.params['num_samples']}
			attr_data = {**attr_data, **attr_features} # merge the two dicts
			self.attractors += [Attractor(self.params, self.G, attr_data)]

		if self.params['debug']:
			total = cp.sum(cp.array([A.size for A in self.attractors]))
			assert(math.isclose(total,1))


		if util.istrue(self.params,'use_phenos'):
			self.build_phenos()


	def build_phenos(self):
		for A in self.attractors:
			found = False
			for P in self.phenotypes:
				if P.id == A.phenotype:
					P.size += A.size 	
					P.attractors += [A]
					assert(not found) # duplicate phenotype
					found=True 

			if not found:
				if '|' in A.phenotype:
					parts = A.phenotype.split('|')
					inpts, outpts = parts[0],parts[1]
				else:
					inpts, outpts = None, None
				self.phenotypes += [Phenotype(A.phenotype, [A],A.size, inpts, outpts)]


	def update_stats(self, lap_stats):
		# dict_keys(['total_avg', 'slow_var', 'windowed_var', 'total_var', 'avg_std_in_time', 'avg_std_in_time_outputs', 'var_threads', 'windowed_var_input_split', 'windowed_var_input_split_sep', 'input_sep'])
		assert(0) # been too long, need thorough debug
		# original keys: total_var, windowed_var, avg_std_in_time, std_btwn_threads
		if len(self.stats) == 0:
			self.stats = deepcopy(lap_stats)
		else:
			for k in lap_stats.keys():
				self.stats[k] += lap_stats[k]


	def _state_avg_to_trinary(self,x,thresh=.1):
		if x<thresh:
			return 0 
		elif x>1-thresh:
			return 1
		else:
			return 2


	def trinary_io_keys(self, X):
		input_ind = self.G.input_indices()
		output_ind = self.G.output_indices()
		ikey = str([self._state_avg_to_trinary(x) for x in X[input_ind]])
		okey = str([self._state_avg_to_trinary(x) for x in X[output_ind]])
		assert('2' not in ikey) # inputs should not fluctuate
		return ikey, okey

	def build_dominant_pheno(self):
		# returns dict {ikey:dominant_okey}
		count = {}
		for A in self.attractors:
			ikey, okey = self.trinary_io_keys(A.avg)
			if ikey not in count:
				count[ikey] = {}
			if okey not in count[ikey]:
				count[ikey][okey] = A.size
			else:
				count[ikey][okey] += A.size
		dominant = {}
		for ikey in count:
			max_val = 0
			for okey in count[ikey]:
				if count[ikey][okey] > max_val:
					max_val = count[ikey][okey]
					dominant[ikey] = okey

		return dominant 

	def __str__(self):
		s="Attractors of SS=\n"
		s+= str([A.avg for A in self.attractors])
		return s

##################################### One Basin #############################################################

def measure(params, G, population=None, x0=None, attractors=True):
	# overview: runs through a transiet period, then measures attractors in equilibrium
	
	if util.istrue(params,['EA','on']):
		assert(population is not None)

	if x0 is None:
		x0 = get_init_sample(params, G)
	feats = {}
	x_in_attractor = lap.lap(params, x0, G, feats, 'transient', use_attractors=attractors, population=population)
	x_final = lap.lap(params, x_in_attractor, G, feats, 'equilibrium',  use_attractors=attractors, population=population)
	#print("simulate.py: avg at eq =",cp.mean(x_final[:,G.nodeNums['FOS']]))
	
	if attractors:
		steadyStates = SteadyStates(params, G) 
		steadyStates.build_attractors(feats)

	if params['debug']:
		features.debug(params, feats, population, G)
	
	if attractors:
		return steadyStates, feats, x_final
	else:
		return feats, x_final


################################# MISC ########################################

def get_init_sample(params, G):
	
	p = .5 #prob a given node is off at start

	if util.istrue(params,['multicellular','active']) and util.istrue(params,['multicellular','init_same']):
		# all samples start identical! 
		assert(0) # have not used in a long time, requires debug
		assert('inputs' not in params.keys() or len(params['inputs'])==0) # otherwise need to make sure they are the same too
		# also would need to check lap.py handling of inputs, ect
		row = cp.random.choice(a=[0,1], size=(G.n), p=[p, 1-p]).astype(bool,copy=False)
		x0 = cp.tile(row,params['num_samples']).reshape((params['num_samples'],G.n))

	else:	
		x0 = cp.random.choice(a=[0,1], size=(params['num_samples'],G.n), p=[p, 1-p]).astype(bool,copy=False)

	#if 'init' in params.keys():
	#	for k in params['init']:
	#		node_indx = G.nodeNums[k]
	#		x0[:,node_indx] = params['init'][k]
	if 'inputs' in params.keys():
		input_indices = G.input_indices()
		input_sets = G.get_input_sets()
		i=0
		for input_set in input_sets:
			istate_start, istate_end = params['input_state_indices'][i][0], params['input_state_indices'][i][1]
			x0[istate_start:istate_end , input_indices] = cp.array(input_set)
			i+=1
		assert(i==2**len(params['inputs']))

	apply_setting_to_x0(params, G, x0)

	#print('temp in basin:')
	#x0[0] = cp.array([0,0,1,0,1,0,1])
	return x0


def format_id_str(x):
	return ''.join(map(str,x.astype(int)))

def cupy_to_numpy(params,d1,d2=None):
	# if using cupy, not extracting these from GPU will lead to SIGNIFICANT slow down

	if params['cupy']:
		for k in d1.keys():
			d1[k]=d1[k].get()
		if d2 is not None:
			for k in d2.keys():
				if k != 'input_sep':
					d2[k]=d2[k].get()
				else:
					for k2 in d2['input_sep']:
						for i in range(len(d2['input_sep'][k2])):
							d2['input_sep'][k2][i] = d2['input_sep'][k2][i].get()


def replace_row_duplicates(a,val,axis):
	unique = cp.sort(a,axis=axis)
	duplicates = unique[:,  1:] == unique[:, :-1]
	unique[:, 1:][duplicates] = val 
	return unique


def apply_setting_to_x0(params, G, x0):
	# in case params have changed, for example sequential run, need to apply mutations to x0
	if 'mutations' in params.keys():
		for k in params['mutations']:
			x0[:,G.nodeNums[k]] = params['mutations'][k]
	if 'init' in params.keys():
		for k in params['init']:
			x0[:,G.nodeNums[k]] = params['init'][k]

	return x0 


def recreate_xf(params, G, SS):
	xf = cp.zeros((params['num_samples'], G.n))
	i=0
	for A in SS.attractors:
		num = int(A.size*params['num_samples'])
		if i+num > params['num_samples']:
			print("\nWARNING simulate.py recreate_xf(): seems to be more samples to recreate than the original number...\n")
			break
		xf[i:i+num] = cp.tile(A.state,(num,1))
		i+=num 
	return xf

