


import re, os, sys

def process_folder(input_folder, output_folder):
	for filename in os.listdir(input_folder):
		input_path = os.path.join(input_folder, filename)
		output_path = os.path.join(output_folder, filename)
		if os.path.isfile(input_path):
			print(f"Processing file: {input_path}")
			process_file(input_path, output_path)

def process_file(input_path, output_path):
	defined_vars = set()
	used_vars = set()
	lines = []

	# Pattern: match variable names (alphanumeric + underscores) not part of operators
	identifier_pattern = re.compile(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b')

	with open(input_path, 'r') as file:
		for line in file:
			line = line.strip()
			if not line:
				continue
			lines.append(line)
			lhs, rhs = line.split(' = ', 1)
			defined_vars.add(lhs.strip())

			vars_in_rhs = identifier_pattern.findall(rhs)
			used_vars.update(vars_in_rhs)

	# Identify variables used but not defined
	undefined_vars = sorted(used_vars - defined_vars)
	remove = ['and','or','not']
	for k in remove:
		if k in undefined_vars:
			undefined_vars.remove(k)

	# Write original lines plus new rows for undefined vars
	with open(output_path, 'w') as file:
		for line in lines:
			file.write(line + '\n')
		for var in undefined_vars:
			file.write(f'{var} = {var}\n')


if __name__ == "__main__":
	if len(sys.argv) !=3:
		sys.exit("Usage: python3 add_inputs.py INPUT_FOLDER OUTPUT_FOLDER")
	process_folder(sys.argv[1], sys.argv[2])