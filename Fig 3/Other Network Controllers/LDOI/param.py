'''
parses parameter file 						
last update 7/17/23	

'''

import os, sys, yaml, util, math, itertools

CUPY, cp = util.import_cp_or_np(try_cupy=0) #should import numpy as cp if cupy not installed

def load(param_file):
	check_file(param_file,'parameter')
	with open(param_file,'r') as f:
		params = yaml.load(f,Loader=yaml.SafeLoader)

	clean(params)
	if 'setting_file' in params.keys():
		assert(0) # depreciating..
		params= load_setting_file(params) # apparently reassigning params within file does not update unless explicitly returned
	elif 'inputs' in params.keys() or 'outputs' in params.keys():
		adjust_for_implicit_setting(params)
	else:
		params['input_state_indices'] = [(0,params['num_samples'])]

	return params

def check_file(file_path,name):
	if not os.path.isfile(file_path):
		sys.exit("Can't find " + name + " file: " + file_path)
	if os.path.splitext(file_path)[-1].lower() != '.yaml':
		sys.exit(name + " file must be yaml format")

def clean(params):
	for k in params.keys():
		param_pow(params, k)
	if 'PBN' in params.keys():
		for k in params['PBN'].keys():
			param_pow(params['PBN'],k)

	CUPY, cp = util.import_cp_or_np(try_cupy=1) #test import
	#params['cupy'] = CUPY # issue: this is always false due to dumb CUPY import

	if util.istrue(params,['PBN','float_update']):
		assert(params['PBN']['active']) # float_update requires PBN
		assert(params['update_rule']=='Gasync') # float_update only makes sense for Gasync


def param_pow(params,k):
	# yaml doesn't maintain json's 10e5 syntax, so here is support for scientific notation. Syntax: 10^5
	if isinstance(params[k], str) and '^' in params[k]:
		parts = params[k].split('^')
		params[k] = int(parts[0])**int(parts[1])

def load_setting_file(params):
	if not os.path.isfile(params['setting_file']):
		sys.exit("Can't find model_file: " + params['setting_file'] + ', check path in parameter file.')
	if os.path.splitext(params['setting_file'])[-1].lower() != '.yaml':
		sys.exit("'setting_file' must be yaml format.")
	
	with open(params['setting_file'],'r') as f:
		model = yaml.load(f,Loader=yaml.SafeLoader)


	if params['debug']:
		shared_items = {k: params[k] for k in params if k in model}
		assert(len(shared_items)==0) #should not have overlaping keys between params and model files

	params = {**model, **params} # in python3.9 can use params_model | params, but may require some people to update python

	adjust_output_thesholds(params)
	reset_settings_for_inputs(params)
	return params

def adjust_output_thesholds(params):
	# can just set single output threshold if they are all the same
	if 'outputs' in params and 'output_thresholds' in params and len(params['outputs'])> 1 and len(params['output_thresholds'])==1:
		params['output_thresholds'] = [params['output_thresholds'][0] for _ in range(len(params['outputs']))]
	elif 'outputs' in params and 'output_thresholds' not in params:
		params['output_thresholds'] = [.5 for _ in range(len(params['outputs']))]


def adjust_for_implicit_setting(params):
	adjust_output_thesholds(params)
	reset_settings_for_inputs(params)
	return params

def reset_settings_for_inputs(params):
	adjust_for_inputs(params)
	add_input_state_indices(params)

def adjust_for_inputs(params):
	if 'inputs' in params.keys():
		k = len(params['inputs'])
		actual_num_parallel = round(params['num_samples']/(2**k))*2**(k)

		if params['num_samples'] < 2**k:
			sys.exit("\nERROR: num_samples parameter must be >= # input combinations, since inputs are run in parallel!\n")

		if actual_num_parallel!=params['num_samples']:
			#if params['verbose']:
			#	print("\nWARNING params.py: num_samples set to",actual_num_parallel,' for equal number of initial states in each input state.')
			params['num_samples']=actual_num_parallel


def add_input_state_indices(params):
	# USED IN massSim.py
	# ie which threads handle which input states
	# not to be confused with the node indices for each input node within a thread
	input_sets = list(itertools.product([0,1],repeat=len(params['inputs'])))
	i=0
	state_indices = []
	for input_set in input_sets:
		state_indices += [(int(i*params['num_samples']/(2**len(params['inputs']))),int((i+1)*params['num_samples']/(2**len(params['inputs']))))]
		i+=1
		
	assert(state_indices[-1][1]==params['num_samples'])
	params['input_state_indices'] = state_indices

def default_output_thresholds(params):
	params['output_thresholds'] = [.5 for _ in range(len(params['outputs']))]

def blank_settings(params, adjust_for_inputs = True):
	# if load one setting, but want to reset to a default
	# note that params['inputs'] and params['outputs'] still need to be correctly re-specified
	default_output_thresholds(params)
	params['init'] = {}
	params['mutations'] = {}
	params['pheno_color_map'] = {} 

	if adjust_for_inputs:
		reset_settings_for_inputs(params)


