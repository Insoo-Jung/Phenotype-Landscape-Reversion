'''
finds controllers using LDOI
implementation of LDOI leverages gpu

'''

INPUT_FOLDER = './models/bnet_exp/'
MUTATION_FILE = './output/38nets_mutation_info.csv'
NOMINAL_ATTRACTOR_FILE = './output/100xT_withoscils_10000S.txt'

import sys, itertools, csv, os, ast
from timeit import default_timer as timer
import util, ldoi_aug, net, major_attractors, param
CUPY, cp = util.import_cp_or_np(try_cupy=1) #should import numpy as cp if cupy not installed

def process_folder(output_file):


	tstart = timer()
	assert(os.path.exists(os.path.dirname(output_file))) # checking that output file can be written
	with open(output_file, 'w'):
		pass  # init file
	params = param.load('./params.yaml') 
	mutations = extract_mutation_data(MUTATION_FILE)
	nominal_attractors = extract_major_attractor_data(NOMINAL_ATTRACTOR_FILE)

	for k in mutations.keys():
		assert(k in nominal_attractors.keys())

	for filename in os.listdir(INPUT_FOLDER):
		input_path = os.path.join(INPUT_FOLDER, filename)
		if os.path.isfile(input_path) and filename in mutations.keys(): # networks without input and output were removed from mutations
			print(f"Processing network: {filename}")
			mutation = mutations[filename]
			nominal_attractor = cp.array(nominal_attractors[filename]['attractor_avg'])
			process_network(input_path, output_file, params, filename, mutation, nominal_attractor)

	tend = timer()
	print("\nTime elapsed: ", round((tend-tstart)/3600,3),'hours')


def process_network(input_path, output_file, params, filename, mutation, nominal_attractor):

	# expand nominal attractor to 2n
	nominal_attractor_expd = cp.concatenate((nominal_attractor,1-nominal_attractor)) # add the not nodes
	get_default_params(params)
	G = net.ParityNet(params, model_file=input_path,read_file_directly=True)
	#val = int(1-nominal_attractor[mutation[0]])
	params['mutations'][G.nodeNames[mutation[0]]] = int(mutation[1])
	G.prepare(params)
	L = ldoi_aug.find_driven(G, params)
	scores = cp.sum((L==1) & (nominal_attractor_expd[None]==1),axis=1)
	controller = int(cp.argmax(scores))
	if controller > G.n: 
		control_val = 0 
		controller = controller - G.n 
	else:
		control_val = 1
	with open(output_file, 'a') as f:
		f.write(filename.replace('_converted.bnet','') + ',\t' + str((controller, control_val)) + '\n')

def get_default_params(params):
	default_params = {
		# generally don't change these
		'time-reversal' : 0,
		'debug':1,
		'init': {},
		'outputs':[],
		'clause_bin_size':1,
		'num_samples':1,
		'mutations':{},
		'init':{},
		'inputs':[],
		'use_inputs':0,	
		'verbose_ldoi':1,
	}
	for k in default_params.keys():
		params[k] = default_params[k]

def extract_mutation_data(mutation_file):
	mutations = {}
	with open(mutation_file, mode='r', newline='', encoding='utf-8') as file:
		reader = csv.DictReader(file)  
		for row in reader:
			network = row['Network_name'].replace('.txt','_converted.bnet')
			mutations[network] = (int(row['Mut_index'])-1, int(row['Mut_Value'])) # -1 converts from R to py indexing
	return mutations

def extract_major_attractor_data(nominal_attractor_file):
	nominal_attractors = {}
	with open(nominal_attractor_file, 'r', encoding='utf-8') as file:
		for line in file:
			line = line.strip()
			if not line:
				continue  
			network, data = line.split(',\t', 1)
			nominal_attractors[network] = ast.literal_eval(data)
	return nominal_attractors


if __name__ == "__main__":
	if len(sys.argv) !=2:
		sys.exit("Usage: python3 ldoi_comparison.py OUTPUT_FILE")
	process_folder(sys.argv[1])