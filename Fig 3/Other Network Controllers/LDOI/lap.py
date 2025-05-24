'''
lap is called by simulate.py to update the network for many time steps						
last update 7/17/23	


'''

import util, sys
import features
CUPY, cp = util.import_cp_or_np(try_cupy=1) #should import numpy as cp if cupy not installed


def lap(params, x0, G, feats, when, use_attractors=False, population=None):

	assert(when in ['transient','equilibrium'])
	x, x0, node_dtype, fixed_nodes = init_arrays(params, G, x0, when) 
	xprev = init_features_for_lap(params, G, feats, population, x0, when, node_dtype, use_attractors)

	for i in range(params['time_steps']):
		#if when=='equilibrium' and i>995 and 0:
		#	print("lap:",x[:,G.nodeNums['Proliferation']].astype(int))
		#	print('\nlap:\n',x[:8].astype(int))
		#print("\n\n\nbefore step:",cp.mean(x,axis=0))
		x = step(params, x, G, fixed_nodes)
		
		features.update(feats, x, xprev, i, when) # feats may need which_nodes at some point, which is curr in step()

		xprev = cp.roll(xprev,-1,axis=0)
		xprev[-1] = x.copy()

	features.finish(params, feats, when)
	return x


def step(params, x, G, fixed_nodes):
	# see junk.py oldlapstep() for other numpy implementations

	node_dtype = util.get_node_dtype(params)
	nodes_to_clauses = G.Fmapd['nodes_to_clauses']        
	clauses_to_threads = G.Fmapd['clauses_to_threads'] 
	threads_to_nodes = G.Fmapd['threads_to_nodes']

	if util.istrue(params,['multicellular','active']):
		# gets expected value of external nodes based on neighbors
		# then sets ON with probability proportional to the expected value
		x_external=(cp.roll(x.astype(cp.float16),1,axis=0)+cp.roll(x.astype(cp.float16),-1,axis=0)+x.astype(cp.float16))/3
		# note that float16 is used since with 3 nghs, each 0 or 1, just need 1/3, 2/3
		
		rd_external = cp.random.random(x.shape)
		x_external_on = rd_external < x_external
		x = (~G.external[cp.newaxis,:] & x) | (G.external[cp.newaxis,:] & x_external_on)
	

	X = cp.concatenate((x,cp.logical_not(x[:,:G.n])),axis=1) # add the not nodes
	clauses = cp.all(X[:,nodes_to_clauses],axis=2) #this is gonna to a bitch to change when clause_index has mult, diff copies

	#print("clauses =\n",cp.mean(clauses,axis=0))
	# partial truths are sufficient for a node to be on (ie they are some but not all clauses)
	partial_truths = cp.any(clauses[:,clauses_to_threads],axis=3)
	#print("partruths =\n",cp.mean(partial_truths,axis=0))

	x_next = cp.sum(cp.matmul(cp.swapaxes(partial_truths,0,1),threads_to_nodes),axis=0).astype(node_dtype)		


	x, which_nodes = apply_update_rule(params, G, x, x_next)
	return x



######################################################################################################

def init_features_for_lap(params, G, feats, population, x0, when, node_dtype, use_attractors):
	features.init(params, G, feats, population, x0, when) 
	twindow = features.get_time_window(feats, when)
	xprev= cp.zeros((twindow,) + x0.shape,dtype=node_dtype)

	features.add_default_features(params, G, feats, x0, when) # incld avg and attractor_ids

	return xprev

def apply_update_rule(params,G, x, x_next):
	if params['update_rule']=='sync':
		x = x_next
		which_nodes = cp.ones(x.shape)
	elif params['update_rule']=='psync': # generalized async, note that this is NOT the same as general async
		p=.5
		which_nodes = cp.random.rand(params['num_samples'], G.n) > p
		x = (which_nodes & x_next) | ((~which_nodes) & x) 
	elif params['update_rule']=='async':
		# WARNING: this is incredibly inefficient (uses an array to update a single element)
		# only use to check against psync		
		which_nodes = cp.zeros((params['num_samples'], G.n), dtype=bool)
		indx1 = cp.arange(params['num_samples'])
		indx2 = cp.random.randint(0, high=G.n, size=params['num_samples'])
		which_nodes[indx1,indx2]=1
		x = (which_nodes & x_next) | ((~which_nodes) & x)
	else:
		sys.exit("\nERROR 'update_rule' parameter not recognized!\n")
	return x, which_nodes

def init_arrays(params, G, x0, when):
	node_dtype = util.get_node_dtype(params)
	x = cp.array(x0,dtype=node_dtype).copy()
	x0 = cp.array(x0,dtype=node_dtype).copy()
	fixed_nodes = None
	return x, x0, node_dtype, fixed_nodes

def debug_print_outputs(params,G,x,target=None):
	k = len(params['inputs'])
	inpt = [G.nodeNums[params['inputs'][i]] for i in range(k)]
	k2 = len(params['outputs'])
	outpt = [G.nodeNums[params['outputs'][i]] for i in range(k2)]
	if target is not None:
		for row in x:
			if cp.all(cp.logical_not(row[inpt]-cp.array(target))):
				print(target,'->',row[outpt].astype(int))
	else: 
		for row in x:
			print(row[inpt].astype(int),'->',row[outpt].astype(int))
