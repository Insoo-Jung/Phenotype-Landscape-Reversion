

import re, os, sys

def process_folder(input_folder, output_folder):
	for filename in os.listdir(input_folder):
		input_path = os.path.join(input_folder, filename)
		output_path = os.path.join(output_folder, filename.replace('.txt','.bnet'))
		if os.path.isfile(input_path):
			print(f"Processing file: {input_path}")
			process_file(input_path, output_path)

def process_file(input_path, output_path):
	lines = []
	with open(input_path, 'r') as file:
		for line in file:
			line = line.strip()
			if not line:
				continue
			line = line.replace(' = ',',\t')
			for pair in [('and','&'),('or','|'),('not','~')]:
				for prior in [' ','(',')','\t']:
					for after in [' ','(',')']:
						line = line.replace(prior + pair[0] + after,prior + pair[1] + after)
			lines.append(line)
	with open(output_path, 'w') as file:
		for line in lines:
			file.write(line + '\n')


if __name__ == "__main__":
	if len(sys.argv) !=3:
		sys.exit("Usage: python3 rewrite_to_bnet.py INPUT_FOLDER OUTPUT_FOLDER")
	process_folder(sys.argv[1], sys.argv[2])