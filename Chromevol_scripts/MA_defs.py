# IMPORTS
import re
import os
import pandas as pd
import numpy as np
import scipy.stats as sc
import sys,argparse,platform
from ete3 import Tree

def fixed_vars():
	CE_res_filename = "/chromEvol.res"  # this name needs to be concatenated to the model's name as a directory
	expectation_file = "/expectations.txt"
	mlAncTree = "/mlAncestors.tree"
	root_freq_filename = "/root_freq"
	sim_control = "/param_sim"
	statistics_names = ["Variance","Entropy","Range","Unique_counts","Parsimony"]
	return(CE_res_filename,expectation_file,mlAncTree,root_freq_filename,sim_control,statistics_names)

### ARGS
def get_arguments():
	parser = argparse.ArgumentParser(description='tests model adequacy of selected ChromEvol model')
	parser.add_argument('--main_res_dir', '-c', help='home directory of ChromEvol results',required=True)
	parser.add_argument('--model_name', '-m', help='the model tested for adequacy', required=True)
	parser.add_argument('--job_id', '-id', help='job ID or Genus name', required=True)
	parser.add_argument('--num_of_trees', '-nt', help='number of trees in the original ChromEvol run',required=False, default=1)
	parser.add_argument('--sims_per_tree', '-ns', help='number of simulations for the adequacy test per original tree',required=False,default=1)
	parser.add_argument('--CE_path', '-ce', help='ChromEvol executable path',required=False, default = "/groups/itay_mayrose/itaymay/code/chromEvol/chromEvol_source-current/chromEvol")
	parser.add_argument('--params', '-p', help='parameters file from user for simulations', required=False, default="")
	parser.add_argument('--counts', '-co', help='counts file', required=True)
	parser.add_argument('--tree', '-t', help='counts file', required=True)
	parser.add_argument('--sanity', '-s', help='counts file', required=True)
	parser.add_argument('--results', '-r', help='counts file', required=True)

	# parse arga
	args = parser.parse_args()
	id = args.job_id
	main_res_dir = args.main_res_dir # + "/" + id + "/"
	in_model = args.model_name
	in_model = in_model.split(",")
	num_of_trees = int(args.num_of_trees)
	sims_per_tree = int(args.sims_per_tree)
	CE_path = args.CE_path
	params_from_user = args.params
	counts_file = args.counts
	tree_full_path = args.tree
	sanity_flag = int(args.sanity)
	results_flag = int(args.results)
	return(id,main_res_dir,in_model,num_of_trees,sims_per_tree,CE_path,params_from_user,counts_file,tree_full_path, sanity_flag,results_flag)