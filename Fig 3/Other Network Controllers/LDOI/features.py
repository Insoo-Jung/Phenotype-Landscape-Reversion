
import util
from itertools import product, combinations
from copy import deepcopy
CUPY, cp = util.import_cp_or_np(try_cupy=1) #should import numpy as cp if cupy not installed


# each possible id specifies a start, update, periodic_update, and finish function	
# avg is mandatory, not sure how will implmt
# jp map population and params directly to each feature?
# TODO: remember to call FeatureSet.debug() from lap.py or smthg


######################## featureS DICT #########################

# functions for features (dict), used to be under FeatureSet class

def init(params, G, feats, population, x0, when):
	last_used = None
	if 'features' not in params:
		return []
	for i in range(len(params['features'])):
		if params['features'][i]['active'] and params['features'][i]['when']==when:
			feats[params['features'][i]['id']] = Feature(params['features'][i], feats, G, params, population)
			if last_used is not None and params['debug']:
				assert(feats[params['features'][i]['id']].features == feats[params['features'][last_used]['id']].features) 
			last_used = i
	start(feats,x0,when)

	if params['debug']:
		for k in feats.keys():
			assert(params is feats[k].params)

	return feats 

def print_features(feats):
	for feature in feats.values():
		if feature.args['print']:
			print(feature)

def add_default_features(params, G, feats, x, when):
	if 'features' not in params.keys():
		params['features'] = []
	for feature_id in ['avg', 'attractor_id']:
		found = False
		for f in params['features']: # don't want to overwrite if already included
			if f['id'] == feature_id:
				found = True 
		if not found:
			params['features'] += [{'id':feature_id ,'active':True, 'print':False,'when':'equilibrium', 'attractor': True, 'update_freq':1}] 
			feats[params['features'][-1]['id']] = Feature(params['features'][-1],feats, G, params, None) # note population not needed for attractor_id and avg

	start(feats, x, when)

def start(feats, x, when):
	for feature in feats.values():
		if feature.when == when:
			feature._start(x)
	
def update(feats, x, xprev, lap_step, when):
	for feature in feats.values():
		if feature.when == when:
			feature._update(x, xprev, lap_step)
	
def finish(params, feats, when):
	# finalizes non-derived features first
	for feature in feats.values():
		if feature.when == when:
			if not util.istrue(feature.args,'derived_feature'):
				feature._finish()
	for feature in feats.values():
		if feature.when == when:
			if util.istrue(feature.args,'derived_feature'):
				feature._finish()

def get_time_window(feats, when):
	# get number of previous time steps that must be kept in memory
	# when can be 'transient' or 'equilibrium'
	time_window = 1
	for feature in feats.values():
		if feature.args['when'] == when and 'update_freq' in feature.args:
			if 'time_window' in feature.args:
				time_window = max(time_window, feature.args['time_window'])
	return time_window

def debug(params, feats, population, G):
	# checking that all point to the same params and population
	# in case, for example, params is modified (NOT reassigned)
	prev=None
	for feature in feats.values():
		assert(feature.params is params)
		assert(feature.population is population)
		assert(feature.G is G)
		if prev is not None:
			assert(feature.params is prev.params)
			assert(feature.population is prev.population)
			assert(feature.G is prev.G)
		prev = feature


############################# Feature Class ######################################

class Feature:
	def __init__(self, args, feats, G, params, population):
		# features is just a dict containing all features
		self.features = feats
		self.args = args
		self.id = args['id']
		self.when = args['when']
		self.val = None

		self.G = G
		self.params = params
		self.population = population

	def __str__(self):
		if self.val.dtype==bool:
			val = self.val.astype(int)
		else:
			val = self.val
		s='\nfeature ' + str(self.id) + ' with shape ' + str(self.val.shape) + ' = \n' + str(val)
		return s

	def _start(self, x):
		# must specify self.val 
		# where x is a 2D matrix of states from lap.py or simulate.py
		if self.id == 'avg':
			self.val = cp.zeros(x.shape) 
		elif self.id in ['oscil','oscil_nodes','oscil_attrs']:
			self.val = cp.zeros(x.shape, dtype=bool) 
		elif self.id in ['flips']:
			shape = (self.args['number_of_flips'],) + x.shape
			self.val = cp.zeros(shape,dtype=bool)
		elif self.id == 'avg_periodish':
			self.val = cp.zeros(x.shape)
			self.last_seen = cp.zeros(x.shape)
			#self.origin = deepcopy(x)  # issue! the x passed here is shape, but not actual starting state!
		elif self.id in ['period', 'no_period']:
			self.val = cp.zeros((self.params['num_samples'],self.params['time_steps']))
			self.last_seen = cp.zeros(self.params['num_samples'])

		elif self.id == 'avg_freq':
			self.val = cp.zeros(x.shape)
			self.last_seen = cp.zeros(x.shape)

		elif self.id == 'px': 
			self.val = cp.zeros(self.size) 
			assert(len(self.val)==len(self.population.genotype))
		elif self.id == 'py': 
			self.val = cp.zeros(len(self.population.targets))
		elif self.id == 'pxy': 
			self.val = cp.zeros((self.population.size,len(self.population.targets)))

		elif self.id == 'attractor_id':
			self.val = x.copy()
			assert(self.val.dtype==bool)
			self.index_dtype = util.get_uint_dtype(self.G.n)
			self.thread_dtype = util.get_uint_dtype(self.params['num_samples'])

			# not sure if pre declaring these helps at all
			# a cleaner broadcast in update() would be more helpful
			self.larger = cp.zeros((self.params['num_samples']),dtype=cp.uint8)
			self.diff = cp.zeros((self.params['num_samples'],self.G.n),dtype=cp.uint8)
			self.first_diff_col = cp.zeros((self.params['num_samples']),dtype=self.index_dtype)
		

		elif self.id in ['temporal_px', 'temporal_pxy']: # probability in time
			# although actually uses 1/2 time window for p_past, p_future
			# assumes ergodic, ie p_past = p_future
			# measures p_past, p_future, p_past&future
			assert(x.dtype==bool) # ie no PBN w floats or probability doesn't work
			if self.id == 'temporal_pxy':
				self.time_window = self.args['time_window']
			elif self.id == 'temporal_px':
				self.time_window = int(self.args['time_window']/2) 
				assert(self.args['time_window'] % 2 == 0) # time window must be even, since also uses 1/2 time iwndow
			self.pindex = cp.array(list(product([0,1],repeat=self.time_window)),dtype=cp.int8)
			self.val = cp.zeros((self.G.n,2**int(self.time_window),self.params['num_samples']))

		elif self.id in ['Ipred_over_time']:
			assert(self.args['time_window']%2==0) 
			self.time_window = self.args['time_window']
			self.half_time_window = int(self.args['time_window']/2) 
			self.pindex = cp.array(list(product([0,1],repeat=self.time_window)),dtype=cp.int8)
			self.half_pindex = cp.array(list(product([0,1],repeat=self.half_time_window)),dtype=cp.int8)

			self.val = cp.zeros((self.G.n,self.params['time_steps']-self.time_window))
			self.H = cp.zeros((self.G.n,self.params['time_steps']-self.time_window))


		### derived features ###
		elif self.id == 'Ipred':
			assert('temporal_px' in self.features.keys() and 'temporal_pxy' in self.features.keys())

	def _update(self, x, xprev, lap_step):
		# note that xprev needs to be indexed sT xprev[-1] is most recent time step
		if 'update_freq' in self.args and lap_step % self.args['update_freq'] == 0: 
			if self.id == 'avg':
				self.val += x
			elif self.id == 'px': 
				assert(x.dtype==bool) # in general should have, also ~x won't work properly otherwise (ie no PBN w floats ect)
				x=cp.append(x,~x,axis=1)
				assert(len(x[0])==len(self.population.genotype[0])) #both should be 2n ie incld expanded nodes
				# would be reasonable to set len(x)=len(genotype) ie same population sizes, but not enforced yet

				# approach:
				#	broadcast genotype & x to match, then let x[1-genotype]=1 such that ignore entries not indexed by genotype 
				#	then px requires that all remaining xs in a row are on (joint probability)
				#	then sum along instances of x

				n,k,m,s = len(x[0]),len(x),len(self.population.targets),len(self.population.genotype)

				# TODO: change tile to cleaner broadcast operation (will save lotto mem)
				x_cast = cp.broadcast_to(x,(s,k,n))
				geno_cast = cp.tile(self.population.genotype,k).reshape(s,k,n)
				# equivalently: geno_cast = cp.transpose(cp.broadcast_to(self.genotype,(k,s,n)),axes=(1,0,2))

				x_masked = x_cast | ~geno_cast
				px_samples = cp.all(x_masked,axis=2)
				self.val +=  cp.sum(px_samples,axis=1)

			elif self.id in ['oscil','oscil_nodes','oscil_attrs']:
				if lap_step > 0:
					self.val = (xprev[-1] != x) | self.val

			elif self.id == 'py': 
				# note that py has to be cast differently than px, since targets are int-indexed
				x=cp.append(x,~x,axis=1)
				py_samples = x[:,self.population.targets]
				self.val += cp.sum(py_samples,axis=0)

			elif self.id == 'flips':
				if lap_step > 0:
					for j in range(1,self.args['number_of_flips']):
						i = self.args['number_of_flips'] - j # need to go in reverse order
						self.val[i] = (self.val[i-1] & (xprev[-1] != x)) | self.val[i]
					self.val[0] = (xprev[-1] != x) | self.val[0]

			elif self.id == 'avg_periodish':
				if lap_step == 0:
					self.origin = deepcopy(x)
				else:
					self.last_seen[self.origin==x] = 0
					self.last_seen[self.origin!=x] += 1
					self.val += self.last_seen

			elif self.id in ['period', 'no_period']:
				if lap_step == 0:
					self.origin = deepcopy(x)
				else:
					returned_to_origin = cp.all(self.origin==x,axis=1)
					self.val[returned_to_origin, lap_step] = self.last_seen[returned_to_origin]+1 # ie the non-zero values are the periods
					self.last_seen[returned_to_origin] = 0
					self.last_seen[~returned_to_origin] += 1

			elif self.id == 'avg_freq':
				if lap_step == 0:
					self.origin = deepcopy(x)
				else:
					self.last_seen[self.origin==x] = 1
					self.last_seen[self.origin!=x] += 1
					self.val += 1/self.last_seen


			elif self.id == 'attractor_id':
				# replace ids with current states that are "larger", where larger is defined by int representation (i.e. first bit that is different)
				self.diff = self.val-x.astype(cp.int8)
				self.first_diff_col = ((self.diff)!=0).argmax(axis=1).astype(self.index_dtype) #i'm worred that this first goes to 64 then converts to 8..
				w = cp.tile(self.first_diff_col,(self.params['num_samples'],1)).T # poss faster if broadcast
				self.larger = cp.take_along_axis(self.diff,w,1)[:,0] # was numpy, change back to np is error
				# gotta be a more efficient way to do this shit, jp need to redo w broadcast
				self.val = (self.larger==-1)[:,cp.newaxis]*x + (self.larger!=-1)[:,cp.newaxis]*self.val # if x is 'larger' than current val, then replace val with it 

			# oscil stuff:
			#if lapStep > 0:
				# isOscil = ~cp.all(xprev==x, axis=1) | isOscil
				# isOscil_indiv = (xprev!=x) | isOscil_indiv

			elif self.id in ['temporal_px','temporal_pxy']:
				assert(xprev.dtype==bool) # can rm soon
				# xprev is a 3D obj shape (t_window, num_samples, num_nodes)

				if lap_step > self.time_window:		
					self.val += cp.all(xprev.T[:,cp.newaxis,:,-self.time_window:]==self.pindex[cp.newaxis,:,cp.newaxis,:],axis=3) # resulting shape: (num_nodes, 2^twindow, num_samples)
					#self.val += cp.sum(p_check, axis=2)
					# note that sliding window occurs anyway due to checking every time step, so just taking first half
			
			elif self.id == 'Ipred_over_time':

				if lap_step >= self.time_window:		
					# need time_window, half_time_window, pindex, half_pindex
					# self.val shape = (num_nodes, steps-self.time_window)
					assert(xprev.T[:,cp.newaxis,:,-self.time_window:-self.half_time_window].shape==xprev.T[:,cp.newaxis,:,-self.half_time_window:].shape)
					p_past = cp.sum(cp.all(xprev.T[:,cp.newaxis,:,-self.time_window:-self.half_time_window]==self.half_pindex[cp.newaxis,:,cp.newaxis,:],axis=3), axis=2)/self.params['num_samples'] # resulting shape: (num_nodes, 2^half_twindow)
					p_futr = cp.sum(cp.all(xprev.T[:,cp.newaxis,:,-self.half_time_window:]==self.half_pindex[cp.newaxis,:,cp.newaxis,:],axis=3), axis=2)/self.params['num_samples'] # resulting shape: (num_nodes, 2^half_twindow)
					p_both = cp.sum(cp.all(xprev.T[:,cp.newaxis,:,-self.time_window:]==self.pindex[cp.newaxis,:,cp.newaxis,:],axis=3), axis=2)/self.params['num_samples'] # resulting shape: (num_nodes, 2^twindow)

					ps = [p_past, p_futr, p_both]
					H_past, H_futr, H_both = 0,0,0
					Hs = [0 for i in range(3)]
					for i in range(3):
						ps[i][ps[i]==0]=1 # log2(1) will be 0 anyway
						Hs[i] = cp.sum(-ps[i]*cp.log2(ps[i]),axis=1) 
					self.val[:,lap_step-self.time_window] = Hs[0] + Hs[1] - Hs[2]
					self.H[:,lap_step-self.time_window] = Hs[0]

	def _finish(self):
		assert(self.args['update_freq']==1) # haven't used yet so reqs explicit debug

		if self.id == 'avg':
			self.val /= self.params['time_steps']

			# TODO: proper derived features...
			if 'multicellular_difference' in self.features.keys():
				k=2
				# later generalize k = difference order
				kwise_combos = cp.array(list(combinations(self.val,k)))
				# shape = (samples choose k, k, G.n)

				# later generalize this to kwise, and diff order norms
				pairwise_node_dist = cp.abs(kwise_combos[:,0,:] - kwise_combos[:,1,:])
				# shape = (samples choose k, G.n)
				# if ord != 1, need to change normzn term too!
				pairwise_dist = cp.linalg.norm(pairwise_node_dist,axis=1,ord=1)/self.G.n
				self.features['multicellular_difference'].val = cp.linalg.norm(pairwise_dist,ord=1)/len(pairwise_dist)

		elif self.id in ['px','py','pxy','temporal_px','temporal_pxy']:
			num_steps = self.params['time_steps']
			if 'time_window' in self.args:
				num_steps -= self.time_window+1 # didn't count before building up all of time window
			#self.val /= (self.params['num_samples'] * num_steps) 
			self.val /= num_steps # note this is per sample per node

			assert(cp.all(cp.isclose(cp.sum(self.val,axis=1),1))) # each node's pr should sum to 1
			# note that dividing by samples and steps simultaneously is an ergodic assumption


		elif self.id == 'oscil_nodes':
			self.val = cp.mean(self.val.astype(int),axis=0) # pr oscil per node
		elif self.id == 'oscil_attrs':
			self.val = cp.mean(self.val.astype(int),axis=1) # pr oscil per sample?
		elif self.id == 'flips':
			self.val = cp.mean(self.val.astype(int),axis=1) # pr to flip per node
			assert(self.val.shape==(self.G.n,))
		
		# oscil stuff:
		#avg_oscils = cp.average(isOscil_indiv) # note this is a scalar: what is pr that a random node in a random attractor oscillates?
		#avg_oscils = cp.average(isOscil_indiv,axis=1)
	
		elif self.id == 'avg_periodish':
			self.val = self.val/self.params['time_steps'] # cp.mean(self.val/self.params['time_steps'],axis=0)
			#print('avg_periodish=',self.val)

		elif self.id == 'period':
			avg_periods = cp.zeros(self.params['num_samples'])
			no_period = (cp.all(self.val==0,axis=1))
			if cp.all(self.val==0):
				print("\nWARNING features.py: all periods are undetermined, setting to period feature =",self.params['time_steps']*10,'\n')
			else:
				avg_periods = cp.true_divide(cp.sum(self.val,axis=1), cp.sum(self.val != 0,axis=1))			
			if cp.any(no_period):
				avg_periods[no_period] = self.params['time_steps']*10 # some large number
			self.val = avg_periods

		elif self.id == 'no_period':
			self.val = (cp.all(self.val==0,axis=1))

		elif self.id == 'avg_freq':
			#any_oscil = cp.any(self.val!=1,axis=0)
			self.val = cp.mean(self.val/self.params['time_steps'],axis=0)

		#### derived features ####
		elif self.id == 'Ipred':
			# one per node
			px = self.features['temporal_px'].val.copy()
			pxy = self.features['temporal_pxy'].val.copy()
			px[px==0]=1 # to avoid error on log, note that log1=0 so entropy will ignore anyway
			pxy[pxy==0]=1

			Hxy = cp.sum(-pxy*cp.log2(pxy),axis=1)
			assert(Hxy.shape == (self.G.n,self.params['num_samples']))
			Hx = cp.sum(-px*cp.log2(px),axis=1)
			assert(Hx.shape == (self.G.n,self.params['num_samples']))


			self.val = - Hxy + 2*Hx 
			# note this assume ergodicity! such that Hpast=Hfuture
			# ie that px and pxy measured in equilibrium

			epsilon=.001
			if not (cp.all(self.val >= 0 - epsilon) and cp.all(self.val <= cp.log2(self.features['temporal_pxy'].time_window) + epsilon)):
				#print("\nWARNING: in features.py out of bounds with Ipred:") #\n",self.val)
				#print('max/min =', cp.max(self.val),cp.min(self.val))
				# this can occur if time window is not long enough to capture dynamics
				#assert(0)
				pass
			
			# TODO: clean all this shit
			if 'Ipred_avg' in self.features.keys():
				maxx=1 #cp.log2(self.features['temporal_pxy'].time_window)
				self.features['Ipred_avg'].val = cp.mean(self.val / maxx)
			if 'Ipred_spc' in self.features.keys():
				maxx=cp.log2(self.features['temporal_pxy'].time_window)
				ipred_avg = cp.mean(self.val / maxx)
				oscils_avg = cp.mean(self.val.astype(bool))
				self.features['Ipred_spc'].val = abs(ipred_avg-oscils_avg) 
				#self.features['Ipred_spc'].val = cp.linalg.norm(cp.mean(self.val,axis=1),ord=8)
			# TODO: this is v poor form, but issue is that Ipred needs to be run first
			#	maybe have some sort of 'priority?' sT depth of derived features is arbitrary?
			if 'Ipred_oscils' in self.features.keys():
				# percent of nodes oscil, ie the avg percent of nodes and samples with non-zero Ipred
				self.features['Ipred_oscils'].val = cp.mean(self.val.astype(bool))
			if 'Ipred_node_max' in self.features.keys():
				# node with largest avg Ipred. ie. max ipred, avgd over samples, of any node 
				self.features['Ipred_node_max'].val = cp.max(cp.mean(self.val,axis=1))
			if 'Ipred_sample_max' in self.features.keys(): 
				self.features['Ipred_sample_max'].val = cp.max(cp.mean(self.val,axis=0))
			if 'Ipred_nodes' in self.features.keys(): 
				self.features['Ipred_nodes'].val = cp.mean(self.val,axis=1)
				assert(self.features['Ipred_nodes'].val.shape==(self.G.n,))
			if 'Ipred_oscil_nodes' in self.features.keys():
				# percent of nodes oscil, ie the avg percent of nodes and samples with non-zero Ipred
				self.features['Ipred_oscil_nodes'].val = cp.mean(self.val.astype(bool),axis=1)
				assert(self.features['Ipred_oscil_nodes'].val.shape==(self.G.n,))
			if 'Entropy_temporal_nodes' in self.features.keys():
				self.features['Entropy_temporal_nodes'].val = cp.mean(Hxy,axis=1)
			if 'Entropy_temporal' in self.features.keys():
				self.features['Entropy_temporal'].val = cp.mean(Hxy)
			if 'Entropy_temporal_over_time' in self.features.keys():
				self.features['Entropy_temporal_over_time'].val = cp.mean(Hxy,axis=0)
				print(self.features['Entropy_temporal_over_time'].val.shape, self.G.n)
				assert(0) # this is the #samples, not #time steps
				# prev assumed ergodic for pxy, px...now don't and instead rely on suffic sample size at each time point...

		if 'avg' in self.args.keys():
			if self.args['avg'] == 'nodes':
				self.val = cp.mean(self.val, axis=0)
			else:
				assert(0) # unknown what to avg over

		if self.params['debug']:
			self._check()

	def _check(self):
		if self.id in ['avg','px','py','pxy','temporal_px','temporal_pxy','flips','oscil_nodes', 'oscil_attrs','oscil']:
			check_standard_normalization(self.val)
		if self.id in []:
			check_standard_normalization(self.val)
			assert(self.val.dtype==bool)
		if self.id == 'Ipred_over_time':
			assert(cp.all(self.val >= -.001))
			assert(cp.all(self.val <= self.time_window))


def check_standard_normalization(x):
	# likely mv to util.py
	if x.dtype!=bool:
		if not(cp.all(x > -.0001) and cp.all(x < 1.0001)):
			print('incorrect normzn:',x)
		assert(cp.all(x > -.0001) and cp.all(x < 1.0001))

#############################################################################################


if __name__ == "__main__":
	import sys, simulate, lap
	# just for debugging
	