
import net, util
from copy import deepcopy 
from ldoi_class import Drivers

CUPY, cp = util.import_cp_or_np(try_cupy=1) #should import numpy as cp if cupy not installed

def find_driven(G, params, inputs = None):

	assert(isinstance(G,net.ParityNet))
	if inputs is not None and len(inputs)==0: 
		inputs = None

	base_pinned_nodes = generate_base_pinned_nodes(G, params, inputs = inputs)
	drivers = augmented_ldoi(G, params, base_pinned_nodes, inputs = inputs) 
	#drivers.remove_base_pinned_nodes() # set inputs and mutated nodes to 0
	return drivers.L


##################################################################################

def generate_base_pinned_nodes(G, params, inputs = None):
	base_pinned_nodes = []	
	for name in params['mutations']:
		base_pinned_nodes += [a_base_pinned_node(G, params, name, params['mutations'][name])]
	for name in params['init']:
		base_pinned_nodes += [a_base_pinned_node(G, params, name, params['init'][name])]
	if inputs is not None:
		assert(len(inputs)==len(params['inputs']))
		for i in range(len(inputs)):
			base_pinned_nodes += [a_base_pinned_node(G, params, params['inputs'][i], inputs[i])]

	return base_pinned_nodes

def a_base_pinned_node(G, params, name, val):
	if val==0:
		return G.nodeNums['!'+name]
	else:
		return G.nodeNums[name]


def print_solutions(G, L):
	drivers = convert_solutions(G, L)
	for k in drivers:
		if len(drivers[k]) > 0:
			print(k,' : ',drivers[k])

def convert_solutions(G,ldoi_solns):
	soln_dict = {}
	for i in range(len(ldoi_solns)):
		soln_names = []
		for j in range(len(ldoi_solns[i])):
			if ldoi_solns[i,j]:
					soln_names += [G.nodeNames[j]]

		soln_dict[G.nodeNames[i]] = soln_names

	return soln_dict

##############################################

def augmented_ldoi(G, params, pinned_nodes,inputs=None,verbose=False):

	loop = 0
	drivers = Drivers(G, params, pinned_nodes, inputs=inputs)
	cont = True
	while cont:
		cont = step(drivers,loop,verbose=verbose)
		loop = util.loop_debug(G.n, loop, expo=40) 
	drivers.L = propagate_matrix(drivers, drivers.L, final_propagate=True)
	return drivers


def step(drivers,loop,verbose=False):

	init_step(drivers)
	add_downstream(drivers)
	cont = drivers.check_continue()
	drivers.update_to_next() 
	return cont

def init_step(drivers):
	# WARNING: possible memory leak in cp
	drivers.L_next = drivers.L.copy()

def add_downstream(drivers):

	drivers.L_next = propagate_matrix(drivers, drivers.L_next)

def propagate_matrix(drivers, M, final_propagate=False):
	# only with 'final_propagate' can !X be added to LDOI(X)
	G = drivers.G
	P = ((M|drivers.pins) & (~drivers.pins_compls)).astype(G.intdtype)
	# note P are the nodes that would be pinned if drivers.pins are pinned
	clauses_on = (P@G.nodes2clauses >= G.composite_counts)

	if not final_propagate:
		M = M | ((clauses_on@G.clauses2nodes) & (~drivers.pins_compls))
	else:
		M = M | ((clauses_on@G.clauses2nodes) & (~drivers.block))
	return M
