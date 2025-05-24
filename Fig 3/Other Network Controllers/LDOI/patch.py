
INPUT_FOLDER = './models/bnet_exp'
INPUT_FILE = './output/100xT_withoscils_10000S.txt'
OUTPUT_FILE = './output/100xT_withoscils_10000S_patched.txt'

import sys, itertools, csv, os, ast
from timeit import default_timer as timer
import util, ldoi_aug, net, major_attractors, param

def process_folder():

	assert(os.path.exists(os.path.dirname(OUTPUT_FILE))) # checking that output file can be written
	with open(OUTPUT_FILE, 'w'):
		pass  # init file
	params = param.load('./params.yaml') 
	datas = extract(INPUT_FILE)

	for filename in os.listdir(INPUT_FOLDER):
		input_path = os.path.join(INPUT_FOLDER, filename)
		if os.path.isfile(input_path):
			print(f"Processing network: {filename}")
			G = net.ParityNet(params, model_file=input_path,read_file_directly=True)
			datas[filename]['node_names'] = G.nodeNames[:G.n]
			with open(OUTPUT_FILE, 'a') as f:
				f.write(filename + ',\t' + str(datas[filename]) + '\n')

def extract(input_file):
	datas = {}
	with open(input_file, 'r', encoding='utf-8') as file:
		for line in file:
			line = line.strip()
			if not line:
				continue  
			network, data = line.split(',\t', 1)
			datas[network] = ast.literal_eval(data)
	return datas

if __name__ == "__main__":
	process_folder()