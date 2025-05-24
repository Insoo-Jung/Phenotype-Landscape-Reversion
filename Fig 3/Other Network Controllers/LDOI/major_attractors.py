'''
for each network in the folder:
finds the major fixed point attractor
finds 5 random mutations that change the major attractor

'''

PARAMS = {
	'time_step_multiplier': 100, # time steps = time_step_multiplier * number_of_nodes
	'num_mutations': 0, # number of relevant mutations to find
	'mutated_dist_thresh':1, # how much does mutated attractor have to be different to count as a significant mutation?
	# for example if all nodes are the same in major nominal and major mutated attractor, except one non-mutated node changes from 0 to 1, then the distance is 1
}

import param, net, simulate, util
import sys, os, random
import numpy as np
from timeit import default_timer as timer

def process_folder(input_folder, output_file):
	tstart = timer()
	assert(os.path.exists(os.path.dirname(output_file))) # checking that output file can be written
	with open(output_file, 'w'):
   		pass  # init file
	params = param.load('./params.yaml') 
	for filename in os.listdir(input_folder):
		input_path = os.path.join(input_folder, filename)
		if os.path.isfile(input_path):
			print(f"Processing network: {filename}")
			process_network(input_path, output_file, params, filename)
	tend = timer()
	print("\nTime elapsed: ", round((tend-tstart)/3600,3),'hours')


def process_network(input_path, output_file, params, filename):
	reset_params(params)
	G = net.ParityNet(params, model_file=input_path,read_file_directly=True)
	modify_params(params, G)
	G.prepare(params)
	steadyStates = simulate.measure(params, G)[0]
	data = get_major_attractor(G, steadyStates)
	if PARAMS['num_mutations']>0:
		data['mutations'] = get_relevant_mutations(G, params, data['attractor_avg'], PARAMS['num_mutations'])
	data['attractor_avg'] = list(data['attractor_avg']) # removing numpy for readability
	data['node_names'] = G.nodeNames[:G.n]
	with open(output_file, 'a') as f:
		f.write(filename + ',\t' + str(data) + '\n')

def reset_params(params):
	params['mutations'] = {}
	params['init'] = {}
	params['inputs'] = []
def modify_params(params, G):
	params['time_steps'] = params['num_samples'] = PARAMS['time_step_multiplier']*G.n
	params['mutations'] = {}

def get_major_attractor(G, SSs):
	max_size = 0 
	attractor_state, attractor_avg = None, None
	some_oscil, all_oscil=False, False
	for A in SSs.attractors:
		if A.size > max_size:
			all_oscil = np.all(A.features['oscil']!=0)
			some_oscil = np.any(A.features['oscil']!=0)
			max_size = A.size 
			attractor_state = A.id
			attractor_avg = A.features['avg']
	if all_oscil:
		print("\tWARNING: all nodes in major attractor oscillate")
	return {'attractor_state':attractor_state, 'attractor_avg':attractor_avg, 'basin_size':max_size, 'some_nodes_oscillate':some_oscil, 'all_nodes_oscillate':all_oscil}


def get_relevant_mutations(G, params, nominal_attractor_avg, k):
	loop = 0 
	mutations = {}
	non_oscil_indices = (nominal_attractor_avg == 0) | (nominal_attractor_avg == 1)
	while len(mutations)<k:
		mutation = random.choice(G.nodeNames).replace('!','')
		nodenum = G.nodeNums[mutation]
		val = random.choice([0,1])
		if nominal_attractor_avg[nodenum] != val:
			params['mutations'] = {mutation:val}
			G.prepare(params)
			steadyStates = simulate.measure(params, G)[0]
			mutated_attractor_avg = get_major_attractor(G, steadyStates)['attractor_avg']
			if not non_oscil_indices[nodenum]:
				base_mutation_dist = 0 # mutated node is oscillating in nominal attractor 
			else:
				base_mutation_dist = abs(val-nominal_attractor_avg[nodenum]) # distance due to mutated node itself
			
			# MUTATION DISTANCE is the distance between non-oscillating nodes in the nominal attractor and their average value in the mutated attractor
			# note the nodes that oscillate in the mutated attractor, but not in the nominal attractor, are counted
			mutated_dist =  np.sum(np.abs(mutated_attractor_avg[non_oscil_indices] - nominal_attractor_avg[non_oscil_indices])) - base_mutation_dist # distance in non-mutated nodes
			
			assert(mutated_dist > -.0001)  # should be non-negative outside of rounding errors
			if mutated_dist >= PARAMS['mutated_dist_thresh']:
				mutations[mutation] = val
				#print("debug: mutated vs nominal attr:", mutated_attractor, nominal_attractor)

		if loop == k*10:
			print("\tWARNING: trying mutation #",k*20)
		loop = util.loop_debug(G.n, loop)

	return mutations

if __name__ == "__main__":
	if len(sys.argv) !=3:
		sys.exit("Usage: python3 add_inputs.py INPUT_FOLDER OUTPUT_FILE")
	process_folder(sys.argv[1], sys.argv[2])