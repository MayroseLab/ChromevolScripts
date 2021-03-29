#from Bio import Phylo
import random
import os
import sys
import shutil
import time
import datetime
import re
from ete3 import Tree
from Bio import Phylo
from ChromEvol_defs import *
from collections import OrderedDict
from glob import glob
import sqlite3
import zipfile

#For Model Adequacy:
from MA_utils import *
from MA_defs import *
from data_processing import process_data
from data_processing import best_model
from data_processing import simulations
from analysis import get_stats
from analysis import test_adequacy

CCDB_DATABASE = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/CCDB_1.46.db'
NCBI_DB_DATABASE = '/groups/itay_mayrose/share/ploidb/ploidb_DBs/GenbankDB_grp_Sep19.db'
CHROMEVOL_CODE_PATH = '/bioseq/chromEvol/Chromevol_scripts/'
CHROMEVOL_EXE_PATH = '/bioseq/chromEvol'
DAILY_TESTS_PATH = '/bioseq/bioSequence_scripts_and_constants/daily_tests'

PARAMS_PATH = '/bioseq/chromEvol/power_PARAM_templates'
DEBUG_FLAG=1
BEST_MODEL_TREE_dict = dict()
EMAIL_FALG = "ON"
TREES_NUM = 1

#----------------------------------------------------------------------------------------------------------#
def MA_Test(tree_dir,BestModel,model_adq_flag,param_dict):

#	BestModel, AIC_per_model_dict = summarize_results(working_dir, tree_dir, models_list, param_dict)
	BestModels_perTree_dict[tree_dir] = BestModel
	bestModelsPerTree_f.write("%s," % BestModel)

	# if model_adq_flag == 'OFF':
	#	continue
	# Otherwise perform Model Adequacy:
	f_ModelAd = open(working_dir + '/ModelAd_log.txt', 'a')
	# AdequacyTest
	# This code was taken from Anna R. code for Model Adequacy and edited by Michal
	if model_adq_flag == 'OFF':
		model_adq_list = ['None']
		return 0
	else:
		model = BestModel  # deault ????
		f_ModelAd.write("Model Adequacy param set to: user model selection\n")
		f_ModelAd.write("param_dict['MA_AdequacyTest']: %s\n" % param_dict['MA_AdequacyTest'])
		f_ModelAd.write("BestModel: %s\n" % BestModel)

		model_adq_list = []
		if param_dict['MA_AdequacyTest'] == 'MA_AdequacyTest_user':
			f_ModelAd.write("Model Adequacy param set to: user model selection\n")
			f_ModelAd.write("User model selection is: %s\n" % param_dict['MA_userModelSelect'])
			model_adq_list = [param_dict['MA_userModelSelect']]
		elif param_dict['MA_AdequacyTest'] == 'ChromBestModel':
			f_ModelAd.write("Model Adequacy param set to: chromevol best model\n")
			f_ModelAd.write("Chromevol best model is: %s\n" % BestModel)
			model_adq_list = [BestModel]
		elif param_dict['MA_AdequacyTest'] == 'AllModels':
			f_ModelAd.write("Model Adequacy param set to: chromevol best model\n")
			f_ModelAd.write("Chromevol best model is: %s\n" % BestModel)
			model_adq_list = models_list

		for model in model_adq_list:
			# -c directory of ChromEvol results (where the dirs of the models are, and results_sum is)
			# -m model name
			# -id an identifier (number of genus name, I think it is used only for printing purposes)
			# -ce ChromEvol path
			# -co counts file (usually named XXcounts_edit)
			# -t tree file (usually named tree_1)
			main_res_dir = tree_dir + '/' + model + '/'

			treeNum = re.search("(\d+)$", tree_dir).group(1)

			#model_adeq_py_cmd = 'python /bioseq/chromEvol/Model_Adequecy/main_for_web.py -c ' \
			#					'%s -m %s -id %s ' \
			#					'-nt 1 -ns 10 -ce %s ' \
			#					'-co %s -s 0 -r 0' \
			#					% (main_res_dir, model, JobID, param_dict['_chromevolExe'], param_dict['_dataFile'])
			# START OF NEW MA PIPELINE
			new_tree_file = main_res_dir + "/mlAncestors.tree"
			new_res_file = main_res_dir + "/chromEvol.res"
			model_adeq_py_cmd = 'python /bioseq/chromEvol/Model_Adequacy/main_for_web.py -c ' \
								'%s -t %s -r %s ' \
								'-out %s -n 1000 ' \
								% (param_dict['_dataFile'], new_tree_file, new_res_file, main_res_dir)
			# END OF NEW MA PIPELINE

			print_to_log(model_adeq_py_cmd, Log_f)
			ma_status = os.system(model_adeq_py_cmd)
			print_to_log(str(ma_status), Log_f)

	return model_adq_list

#----------------------------------------------------------------------------------------------------------#
def create_ploidy(working_tree_dirs,working_dir):

	thresFilePP = working_dir + '/thresholds_PP'
	thresFileDP = working_dir + '/thresholds_DP'

	## Determine threshold for calling diploids/polyploids (using simulations)
	ploidyCallType = "DUPL DEMI BASE"  # MD: ???? any other options?

	print_to_log("Before determine_thresholds_for_power line 102\n", working_dir + '/Log.txt')

	thresPP, thresDP = determine_thresholds_for_power(working_tree_dirs, thresFilePP, thresFileDP, ploidyCallType)
	print(thresPP)
	print(thresDP)

	print_to_log("thresPP %f thresDP %f\n"%(thresPP,thresDP), working_dir + '/Log.txt')

	## Compute reliability from real inference
	relInfer = working_dir + '/reliability_infer.txt'
	if compute_reliability_from_real_for_power(working_tree_dirs, thresPP, thresDP, relInfer,
											   ploidyCallType) != 'DONE':
		print("ERROR: Failed to compute reliability from real inferences\n")
	inferReliabilityRef = parse_reliability_file(relInfer)

	## Compute reliability from simulations inference
	relSim = working_dir + '/reliability_sims.txt'
	if compute_reliability_from_sim_for_power(working_tree_dirs, thresDP, thresPP, ploidyCallType, relSim,
											  inferReliabilityRef) != 'DONE':
		print("ERROR: Failed to compute reliability from simulations\n")

	## Summarize reliability -  use counts edit file
	countsEditFile = param_dict[
		'_dataFile']  # MD changed from: my @countsEdit = grep {/counts_edit$/} readdir(DIR);
	finalOut = working_dir + '/ploidy.csv'
	summarize_reliability_for_power(relSim, relInfer, countsEditFile, finalOut)
	print_to_log("End OF create_ploidy functions",working_dir + '/Log.txt')
	return

#----------------------------------------------------------------------------------------------------------#
def send_start_email(param_dict,working_dir):

	if EMAIL_FALG == 'OFF':
		return
	JobID = os.path.basename(os.path.normpath(working_dir))
	email_addr = param_dict['inputEmail']
	jobTitle = param_dict['jobTitle'].replace(" ", "%20")
	os.system('perl %ssendFirstEmail.pl --toEmail %s --id %s --jobTitle "%s"' %(CHROMEVOL_CODE_PATH,email_addr,JobID, jobTitle))
	return
#----------------------------------------------------------------------------------------------------------#
def send_end_email(param_dict,working_dir,status_str):

	if EMAIL_FALG == 'OFF':
		return
	JobID = os.path.basename(os.path.normpath(working_dir))
	email_addr = param_dict['inputEmail']
	jobTitle = param_dict['jobTitle'].replace(" ", "%20")
	if jobTitle != "daily%20test":
		os.system('perl %ssendLastEmail.pl --toEmail %s --id %s --JobStatus %s --jobTitle "%s"' %(CHROMEVOL_CODE_PATH,email_addr,JobID,status_str, jobTitle))
	else:
		date = datetime.datetime.today().strftime('%d%m%Y')
		with open(os.path.join(DAILY_TESTS_PATH, f'chromEvol_{date}.txt'), "w") as f:
			resultsLink = f'http://ChromEvol.tau.ac.il/results.html?jobId={JobID}&jobTitle={jobTitle}'
			f.write(f'{status_str},{resultsLink}')
		f.close()	
		
	return

#----------------------------------------------------------------------------------------------------------#
def copy_dir_to_another_dir(source_dir, target_dir,Log_f):
	path, dirname = os.path.split(source_dir)
	target_dirname = os.path.join(target_dir, dirname)
	print_to_log(("About to copy dir %s to %s" %(source_dir,target_dirname)), Log_f)
	shutil.copytree(source_dir, target_dirname)
#----------------------------------------------------------------------------------------------------------#
def get_counts(counts_path):

	counts_dict = dict()
	with open(counts_path,'r') as counts_f:
		for line in counts_f:
			if '>' in line:
				taxa_name = line.strip().replace('>','')
			else:
				counts_data = line.strip()
				counts_dict[taxa_name] = counts_data
				taxa_name=''
				continue
	return counts_dict
#----------------------------------------------------------------------------------------------------------#
def get_newick_branch_labels(newick_tree_file):

	branch_length_label_dict=dict()
	with open(newick_tree_file,'r') as newick_f:
		for line in newick_f:
			node_length_list = re.findall("(\[N\d+\-\d+\]:\d+\.?\d*)",line)
			break
		for item in node_length_list:
			node_name = item.split(':')[0].replace('[','').replace(']','')
			length_num = item.split(':')[1]
			branch_length_label_dict[length_num]=node_name

	return branch_length_label_dict


#----------------------------------------------------------------------------------------------------------#
def get_ploidy(ploidy_path):

	ploidy_dict = dict()
	with open(ploidy_path,'r') as ploidy_f:
		for line in ploidy_f:
			line_split = line.replace('"','').split(',')
			taxa_name = line_split[0]
			ploidy_inf = line_split[2]
			ploidy_dict[taxa_name]=ploidy_inf

	return ploidy_dict




#---------------------------------------------------------------------------------------------------------------------#
def create_tables_from_expectations(results_dir,exp_file):

	# This function will create the tables for the following:
	#1. #ALL EVENTS EXPECTATIONS PER NODE (ended with #+++++++++++++++++++++++++++++)
	#2. #Expected number of events from root to leaf (at the end of file so blank line is the end)
	blk_events_f = open(results_dir+'/Results_Events.csv','w')
	blk_expNum_f = open(results_dir+'/Results_ExpNum.csv','w')
	BLK_Events_lines = []
	BLK_Num_Events_lines = []
	with open(exp_file,'r') as exp_f:
		copy_blk_events = 0
		copy_blk_Num_events = 0

		for line in exp_f:
			if '#ALL EVENTS EXPECTATIONS PER NODE' in line:
				copy_blk_events = 1
				continue
			if copy_blk_events == 1 and '#+++++++++++++++++++++++++++++' in line:
				copy_blk_events = 0
				continue
			if copy_blk_events == 1:
				blk_events_f.write(line.strip().replace('\t',','))
				blk_events_f.write('\n')
				#BLK_Events_lines.append(line.strip())
			if '#Expected number of events from root to leaf' in line:
				copy_blk_Num_events = 1
				continue
			if copy_blk_Num_events == 1 and not line.strip():
				copy_blk_Num_events = 0
				continue
			if copy_blk_Num_events == 1:
				blk_expNum_f.write(line.strip().replace('\t',','))
				blk_expNum_f.write('\n')
				#BLK_Num_Events_lines.append(line.strip())

	return

#---------------------------------------------------------------------------------------------------------------------#
def convert_to_PhyD3(input_tree,counts_path,ploidy_path):
#ploidy_path - ploidy_path = dir_path + '/ploidy.csv'
#counts_path = dir_path + '/countsFile'

	phyD3_tree = input_tree.split('.')[0] + '_phylo'
	phyD3_tree_edit = input_tree.split('.')[0] + '_phylo_edit'
	node_length_dict = get_newick_branch_labels(input_tree)
	Phylo.convert(input_tree, 'newick', phyD3_tree, 'phyloxml')

	# Create counts and Ploidy dictionaries:
	taxa_counts_dict = get_counts(counts_path)
	print("convert_to_PhyD3 def: ploidy_path value")
	print(ploidy_path)
	if ploidy_path is not 'None':	#WITH Ploidy !!!
		taxa_ploidy_dict = get_ploidy(ploidy_path)
		for item in taxa_ploidy_dict.keys():
			print(item + ':' + taxa_ploidy_dict[item])
			print("Inside convert_to_PhyD3 taxa ploidy dict data")

	First_line = '<?xml version="1.0" encoding="UTF-8"?>\n'
	if ploidy_path is not 'None':	#WITH Ploidy !!!
		lable_section = '<labels>\n\
							<label type="text">\n\
							  <name show="1">Chromosome Count</name>\n\
							  <data tag="chrom"/>\n\
							</label>\n\
						   <label type="color">\n\
							  <name>Color label</name>\n\
							  <data tag="colortag"/>\n\
							</label>\n\
							<label type="text">\n\
							  <name show="2">Ploidy Inference</name>\n\
							  <data tag="ploidy"/>\n\
							</label>\n\
						 </labels>'
	else:
		lable_section = '<labels>\n\
								<label type="text">\n\
								  <name show="1">Chromosome Count</name>\n\
								  <data tag="chrom"/>\n\
								</label>\n\
							 </labels>'

	with open(phyD3_tree, 'r') as phylo_f:
			with open(phyD3_tree_edit, 'w') as phylo_edit_f:
				phylo_edit_f.write(First_line)
				for line in phylo_f:
					# Get the name of taxa if exists:
					if '<name>' in line:
						taxa_name_withCount = re.split('>|<|\n', line)[2]
						taxa_name = taxa_name_withCount.split('-')[0]
						phylo_edit_f.write(line.replace(taxa_name_withCount, taxa_name))
						if taxa_name in taxa_counts_dict.keys():
							chrom_tag = '<chrom>' + taxa_counts_dict[taxa_name] + '</chrom>\n'
						else:
							chrom_tag = '<chrom> - </chrom>\n'
						phylo_edit_f.write(chrom_tag)
						if ploidy_path is not 'None':	#WITH Ploidy !!!
							print("Inside convert_to_PhyD3 line 300")
							if taxa_ploidy_dict[taxa_name] == '0':
								phylo_edit_f.write('<colortag>0xF7104A</colortag>\n')
							else:
								phylo_edit_f.write('<colortag>0x16db37</colortag>\n')
							ploidy_tag = '<ploidy>' + taxa_ploidy_dict[taxa_name] + '</ploidy>\n'
							phylo_edit_f.write(ploidy_tag)
						continue
					if '<branch_length>' in line:
						branch_length = re.split('>|<|\n', line)[2]
						phylo_edit_f.write(line)
						if branch_length in node_length_dict.keys():
							line_added = '<name>%s</name>' % (node_length_dict[branch_length])
							phylo_edit_f.write(line_added)
						continue
					if '<phylogeny' in line:
						phylo_edit_f.write(line)
						phylo_edit_f.write(lable_section)
						continue
					else:
						phylo_edit_f.write(line)
						continue

#----------------------------------------------------------------------------------------------------------#
def return_min_param_Model(selected_models):

	models_list = selected_models.split(',')
	min_param_cnt = 99999
	min_mod = 'xxx'
	for mod in models_list:
		count_line = 0
		param_file = PARAMS_PATH + '/param_' + mod #param_CONST_RATE
		with open(param_file,'r') as par_f:
			for line in par_f:
				if line:
					count_line+=1
		if count_line < min_param_cnt:
			min_mod = mod
	return min_mod


#----------------------------------------------------------------------------------------------------------#
def return_MA_vecData(vec_file):

	#Adeq_stats_names_VecLocation = {'Variance':0 ,'Entropy':1,'Parsimony':4,'Time_parsimony':5}
	vec_data = dict()
	with open(vec_file,'r') as vec_f:
		# commented out by Josef
		#for line in vec_f: 
		#	line_data = line.strip().replace('[','').replace(']','') # [1, 1, 1, 1, 1, 0]
		#	vec_data = line_data.split(',')
		vec_data = vec_f.readlines()

	#Add summary value:
	count_in_range = 0
	for vec_val in vec_data:
		if float(vec_val) > 2.5 and float(vec_val) < 97.5:
			count_in_range+=1
	if count_in_range == 4:
		vec_data.append("Adequate for all") #XX out of 4 is Adequate" or "Adequate for all" or "Adequate for none
	if count_in_range < 4 and count_in_range > 1:
		vec_data.append("%d out of 4 are Adequate") #XX out of 4 is Adequate" or "Adequate for all" or "Adequate for none
	if count_in_range == 1:
		vec_data.append("1 out of 4 is Adequate") #XX out of 4 is Adequate" or "Adequate for all" or "Adequate for none
	if count_in_range == 0:
		vec_data.append("Adequate for none") #XX out of 4 is Adequate" or "Adequate for all" or "Adequate for none

	return vec_data

#----------------------------------------------------------------------------------------------------------#


def changeModelName(str):
	if str == "CONST_RATE_NO_DUPL":
		return "Dys"
	if str == "CONST_RATE":
		return "DysDup"
	if str == "CONST_RATE_DEMI":
		return "DysDupDem*"
	if str == "CONST_RATE_DEMI_EST":
		return "DysDupDem"
	if str == "BASE_NUM":
		return " DysBnum"
	if str == "BASE_NUM_DUPL":
		return " DysBnumDup"


def Create_results_for_web(results_dir,models_list,model_adq_list,AIC_per_model_dict):

	select_model_file = '/Chromevol_selected_models.txt'
	chrom_res = 'chromEvol.res'
	tree_paths = dict()
	tree_sel_models = dict()
	with open(results_dir + select_model_file, 'r') as models_f:
		for line in models_f:
			# Selected model (AIC) for Tree# 1: BASE_NUM_DUPL,
			line_split = line.split(':')
			selected_model = line_split[1].replace(',', '').strip()
			tree_num = line_split[0].split(' ')[-1]
			tree_path = results_dir + '/ChromEvol_Tree_' + str(tree_num) + '/'
			tree_paths[tree_num] = tree_path
			tree_sel_models[tree_num] = selected_model
			print_to_log("results: tree %s: selected model: %s"%(tree_num,selected_model), results_dir + '/Log.txt')

	out_f = open(results_dir + '/RESULT_PRINTOUT.csv', 'w')
	if model_adq_list:
		out_f2 = open(results_dir + '/MA_RESULT_PRINTOUT.csv', 'w')
	model_params_data_dict = dict()
	param_keys = ['BASE_NUMBER_R', 'DUPL', 'GAIN_CONST', 'HALF_DUPL', 'LOSS_CONST', 'BASE_NUMBER', 'Distance from best AIC']
	for tree_num in tree_paths.keys():	#at this point 1 tree only
		for model_name in models_list:
			print_to_log("results: line 286 - model %s"%model_name, results_dir + '/Log.txt')
			chrom_num_dict = dict()
			gain_loss_dict = dict()
			exp_file = tree_paths[tree_num] + model_name + '/expectations.txt'
			create_tables_from_expectations(tree_paths[tree_num] + model_name, exp_file)
			# Read chromEvol.res file:
			with open(tree_paths[tree_num] + model_name + '/' + chrom_res, 'r') as res_f:
				strt_gain_loss = 0
				for line in res_f:
					if 'chromosome' in line:  # min chromosome in data: 7
						line = line.replace('#', '')
						line_splt = line.split(' ')
						if 'allowed' in line:
							chrom_num_type = line_splt[0] + ' chromosome allowed'
						else:
							chrom_num_type = line_splt[0] + ' chromosome in data'
						chrom_num_dict[chrom_num_type] = line_splt[-1].strip()

					for param_id in param_keys:
						#print_to_log("results: line 306 - param_id %s" % param_id, results_dir + '/Log.txt')
						if param_id in line:
							if param_id == 'BASE_NUMBER' and 'BASE_NUMBER_R' not in line:
								gain_loss_dict[param_id] = int(line.strip().split('\t')[1])
								if (gain_loss_dict[param_id] == 0.0):
									gain_loss_dict[param_id] = '0'
								else:
									gain_loss_dict[param_id] = int(str(gain_loss_dict[param_id]))
							else:
								gain_loss_dict[param_id] = format(float(line.strip().split('\t')[1]),'.4f') # format(a, '.2f')format(a, '.2f')
								if (gain_loss_dict[param_id] == 0.0):
									gain_loss_dict[param_id] = '0'
								else:
									gain_loss_dict[param_id] = str(gain_loss_dict[param_id])
			model_params_data_dict[model_name]=gain_loss_dict


	counts_data_f = open(results_dir+'/chrom_counts_data.txt','w')
	print_to_log("results: line 329 model %s - " % (results_dir+'/chrom_counts_data.txt'), results_dir + '/Log.txt')
	#data_range=''
	#allowed_range=''
	#for count_param in chrom_num_dict.keys():
	#	print_to_log("results: line 333 count_param %s - " % count_param,results_dir + '/Log.txt')
#
	#	if count_param == 'min chromosome in data':
	#		data_range += '%s-'%chrom_num_dict[count_param]
	#	elif count_param == 'max chromosome in data':
	#		data_range += '%s' % chrom_num_dict[count_param]
	#	if count_param == 'min chromosome allowed':
	#		allowed_range += '%s-'%chrom_num_dict[count_param]
	#	elif count_param == 'max chromosome allowed':
	#		allowed_range += '%s' % chrom_num_dict[count_param]
	#counts_data_f.write(" Chromosome in data range: %s , Chromosome allowed range: %s"%(data_range,allowed_range))
	counts_data_f.write(" Chromosome in data range: TEMP , Chromosome allowed range: TEMP")
	counts_data_f.close()

	#Sort models list according to AIC values list:
	model_list_sort =  list({k: v for k, v in sorted(AIC_per_model_dict.items(), key=lambda item: item[1])}.keys())
    
    # Josef: fix sort in case of user selected model
	print_to_log("model_adq_list.len = %d" %len(model_adq_list), working_dir + '/Log.txt')
	if len(model_adq_list) == 1: 
		pos = model_list_sort.index(model_adq_list[0])
		print_to_log("pos = %d" %pos, working_dir + '/Log.txt')
		if pos != 0: 
			model_swap = model_list_sort[0]
			model_list_sort[0] = model_list_sort[pos]
			model_list_sort[pos] = model_swap

	#Write models in file according to AIC order:
	with open(working_dir + '/selected_models.txt','w') as f_models:
		comma = ''
		for model_name in model_list_sort:
			f_models.write(comma+convertToNewModelNames(model_name)) # josef
			comma = ','
	f_models.close()

	out_f.write('"Header1",')
	if model_adq_list: out_f2.write('"Header1",')
	idx = 2
	for model_name in model_list_sort:
		out_f.write('"Header%d",' % idx)
		if model_adq_list:
			if model_name in model_adq_list:
				out_f2.write('"Header%d",' % idx)
		idx += 1
	out_f.write('\n')
	if model_adq_list: out_f2.write('\n')

	out_f.write('"",')
	None_str = "-"
	for model_name in model_list_sort:
		# new_model_name = changeModelName(model_name) # Added by Anna
		out_f.write('"%s",' % model_name) #  old model names
		# out_f.write('"%s",' % new_model_name)
	out_f.write('\n')
	for idx in range(0, len(param_keys)):
		out_f.write('"%-25s",' % param_keys[idx])
		for model in model_list_sort:
			if param_keys[idx] in model_params_data_dict[model].keys():
				out_f.write('"%-15s",' % model_params_data_dict[model][param_keys[idx]])
				continue
			elif param_keys[idx] == 'Distance from best AIC':
				out_f.write('"%-15s",' % str(AIC_per_model_dict[model]))
				continue
			else:
				out_f.write('"%-15s",' %None_str)

		out_f.write('\n')

	#Print the model selection results and Model Adequacy:
	Adeq_stats_names_VecLocation = {'Variance':0 ,'Entropy':1,'Parsimony':2,'Parsimony vs. time':3}

	if model_adq_list: out_f2.write('"",')
	dict_ma_for_models=dict()
	if model_adq_list:
		for model_name in model_list_sort:
			if model_name in model_adq_list:
				out_f2.write('"%s",' % model_name)
				#vec_file = working_dir + '/ChromEvol_Tree_1/' + model_name + '/adequacy_test/adequacy_vec'
				#vec_file = working_dir + '/ChromEvol_Tree_1/' + model_name + '/true_percentiles'
				vec_file = working_dir + '/ChromEvol_Tree_1/' + model_name + '/PVs' # Josef
				model_VecData = return_MA_vecData(vec_file)
				dict_ma_for_models[model_name] = model_VecData
		out_f2.write('\n')
		for MA_param in Adeq_stats_names_VecLocation.keys():
			out_f2.write('"%-25s",' % MA_param)
			for model in model_list_sort:
				if model in model_adq_list:
					out_f2.write('"%-15s",' % dict_ma_for_models[model][Adeq_stats_names_VecLocation[MA_param]])
				else:
					c_spc = ' '
					out_f2.write('"%-15s ",' % c_spc)
			out_f2.write('\n')
		#Statistics line: commented out by Josef
		#out_f2.write('"Statistics distributions",')
		#for model in model_list_sort:
		#	if model in model_adq_list:
		#		out_f2.write('"DistLink",')
		#	else:
		#		out_f2.write('"",')
		#out_f2.write('\n')

	return

#----------------------------------------------------------------------------------------------------------#
def read_job_id(f_job_num):

	with open(f_job_num,'r') as job_f:
		for line in job_f:
			job_id = line.strip()
	return job_id
#----------------------------------------------------------------------------------------------------------#
# creates the base number transitions probability vector for simulations
def create_bntpv(expFile, mlAncTree, baseNum,model_name):

	with open(expFile,'r') as expFile_f:
		perNode_EXPECTATIONS_lines=[]
		look_for_end=0
		for line in expFile_f:
			if "#ALL EVENTS EXPECTATIONS PER NODE" in line and look_for_end == 0:
				look_for_end = 1    #stop when you find the next #+++++++++++++++++++++++++++++ line
			if look_for_end == 1 and "#++++++++++++++++++++++++" in line:
				print("Done copy evenets lines !!!")
				break
			if look_for_end == 1:
				if "#ALL EVENTS EXPECTATIONS PER NODE" not in line and "BASE-NUMBER" not in line:
					perNode_EXPECTATIONS_lines.append(line.strip())
	expFile_f.close()

	btnodes=dict()
	if DEBUG_FLAG == 1: print("create_bntpv model is: %s\n" %model_name)
	if model_name == 'BASE_NUM':
		for item in perNode_EXPECTATIONS_lines:
			item_split = item.split()
			node = item_split[0]
			dupl_weight = float(item_split[3])  # Correction Oct 2019 (if Dupl is 1 then auto we set Base num to 1)
			weight = float(item_split[5])
			if weight > 0.1 or dupl_weight > 0.5:
				if dupl_weight > 0.5:
					weight = dupl_weight
				elif weight > 1:
					weight = 1
				btnodes[node] = weight
		if not btnodes:  # if for some reason there are no base transition events
			if DEBUG_FLAG == 1: print("btnodes is empty so vector is %d=1.00" % baseNum)
			return "%d=1.00" % baseNum
	elif model_name == 'BASE_NUM_DUPL':
		for item in perNode_EXPECTATIONS_lines:
			item_split = item.split()
			node = item_split[0]
			weight = float(item_split[5])
			if weight > 0.1:
				if weight > 1:
					weight = 1
				btnodes[node] = weight
		if not btnodes: # if for some reason there are no base transition events
			if DEBUG_FLAG == 1: print("btnodes is empty so vector is %d=1.00"%baseNum )
			return "%d=1.00"%baseNum

	if DEBUG_FLAG == 1:
		print("btnodes[node_name] printout:\n")
		for node_name in btnodes.keys():
			print("Node: %s, weight: %.2f\n" % (node_name,btnodes[node_name]))

	### for each node, find the base transition
	transitions = dict()  # will contain all transitions found in data (value is # of times transition was found)
	tree_test=Tree(mlAncTree, format=1)
	for node in tree_test.traverse("postorder"):
		nodeid = re.search("\[?([^\-]+)-([^\]]+)\]?", node.name)
		nodeName = nodeid.group(1)
		if nodeid.group(2) is not 'x' and nodeid.group(2) is not 'X':
			nodeCount = int(nodeid.group(2))
		else:
			nodeCount = nodeid.group(2)
		# if node has base num transition
		if nodeName in btnodes.keys():
			ancestor = node.up
			ancestorid = re.search("\[?([^\-]+)-([^\]]+)\]?", ancestor.name)
			ancestorName = ancestorid.group(1)
			if ancestorid.group(2) is not 'x' and ancestorid.group(2) is not 'X':
				ancestorCount = int(ancestorid.group(2))
			else:
				ancestorCount = ancestorid.group(2)
			# if both node and parent have counts
			if nodeCount is not 'x' and nodeCount is not 'X' and ancestorCount is not 'x' and ancestorCount is not 'X':
				#Need to add the posibility of probabilities vector !!!
				transition = nodeCount - ancestorCount  # may include gains\losses
				realTrans = round(transition /baseNum) * baseNum
				# add transition to hash
				if realTrans > 0:
					if realTrans in transitions.keys():
						transitions[realTrans] += btnodes[nodeName]
					else:
						transitions[realTrans] = btnodes[nodeName]
			else:
				btnodes[nodeName] = "NA"

	# in case no transitions were found, give transition by base number prob of 1
	if not transitions:
		transitions[baseNum] = 1.00

	### calculate probabilities
	allTransitions = list(transitions.values())
	totalTrans = sum(allTransitions)
	# remove rare transitions (less then 0.05)
	transitions_copy = transitions.copy()
	for key in transitions.keys():
		if (transitions[key]/totalTrans < 0.05):
			del transitions_copy[key]
	#copy the edited dictionary back to the original one:
	transitions = transitions_copy.copy()

	# create final probs hash
	probsHash=dict()
	for key in transitions.keys():
		res=float(transitions[key]/totalTrans)
		finalProb = "%.2f" %res # limited to 2 decimal digits
		probsHash[key] = float(finalProb)

	# make sure probs sum to 1
	allProbs = list(probsHash.values())
	probSum = sum(allProbs)
	diff = 1-probSum
	if diff != 0:
		# complete or substract to 1 on a random prob
		transitions = list(probsHash.keys())
		randTrans = transitions[0]
		probsHash[randTrans] += diff

	# create probs vector
	vector = ""
	for key in probsHash.keys():
		vector+= str(key) + '=' + str(probsHash[key]) + '_'
	vector = vector[:-1]

	return vector

#----------------------------------------------------------------------------------------------------------#
def test_simulations(simNum,allSimDir):#, maxInData, minInData, simulationDir):

	simulation_report_f = open(allSimDir +'/sim_Report.txt','w')
	n = 1
	maxRange = 200  # max range allowed in simulations
	rangeTest=0;overMaxTest=0;ploidyTest=0

	while n <= simNum*2:
		simDir = allSimDir+'/dir_sim_' + str(n) + '/0/'
		simMSA = simDir + '/simCounts.txt'
		simEvents = simDir + '/simEvents.txt'

		# check if maximum count range was reached
		chromNumbers=[] # doesn't include 'x' counts for min/max calc
		chromNumbers_data=[]    #incdlue x data
		try:
			sim_msa_f = open(simMSA,'r')
			for line in sim_msa_f:
				if not re.search('^(^>)',line):
					chromNumbers_data.append(line.strip())
					chrom_num = int(line.strip())
					chromNumbers.append(chrom_num)
			maxChromNum = max(chromNumbers)
			minChromNum = min(chromNumbers)
			if (maxChromNum - minChromNum < maxRange):
				rangeTest = 1
			else:
				sim_msa_f.close()
				n += 1
				continue
		except IOError:
			simulation_report_f.write("Can't open simulation MSA file %s\n" %simMSA)
		sim_msa_f.close()


		# check if sim reached max chrom allowed
		try:
			sim_ev_f = open(simEvents,'r')
			for line in sim_ev_f:
				if re.search('(Total number of transitions to max chromosome: \d+)',line):
					num_trans = int(re.search('Total number of transitions to max chromosome: (\d+)',line).group(1))
					simulation_report_f.write("num_trans:%s\n"%str(num_trans))
					#simulation_report_f.write(str(num_trans))
			if num_trans == 0:
				overMaxTest = 1
				sim_ev_f.close()
			else:
				sim_ev_f.close()
				n += 1
				continue
		except IOError:
			simulation_report_f.write("Can't open simulations events file %s" %simEvents)
		sim_ev_f.close()

		# check if all leafs are di/poly ploids
		try:
			check_lines=False
			sim_ev_f = open(simEvents,'r')
			relevant_lines=[]
			for line in sim_ev_f:
				if check_lines or 'LEAF	GAIN	LOSS	DUPLICATION	DEMI-DUPLICATION' in line:
					check_lines = True
					relevant_lines.append(line.strip())
			leafTotalNum = 0
			polyPloidNum = 0
			for line in relevant_lines:
				if 'LEAF	GAIN	LOSS	DUPLICATION	DEMI-DUPLICATION' in line:
					continue
				eventLine = line.split('\t')
				#eventLine = re.search("[^\t]+\t\d+\t\d+\t(\d+)\t(\d+)(\t(\d+))?",line)
				leafTotalNum+=1
				if len(eventLine) > 5:    # There is BASE-NUMBER col
					simulation_report_f.write(" There is BASE-NUMBER col\n")
					if (eventLine[3] is not '0' or eventLine[4] is not '0' or eventLine[5] is not '0'):
						polyPloidNum+=1
				else:
					simulation_report_f.write(" There is NO BASE-NUMBER col\n")
					if (eventLine[3].strip() is not '0' or eventLine[4].strip() is not '0'):
						polyPloidNum += 1
			if (polyPloidNum < leafTotalNum and polyPloidNum > 0):
				ploidyTest = 1
		except IOError:
			simulation_report_f.write("Can't open simulations events file %s\n" %simEvents)
		sim_ev_f.close()

		seen=dict()
		uniq=0
		for count in chromNumbers_data:
			if count == 'x':
				continue
			if count not in seen.keys():
				uniq+=1
				seen[count]=1
			else:
				seen[count]+=1
		if uniq == 1:
			ploidyTest = 0

		# decide if sim is OK
		if ploidyTest == 1 and overMaxTest == 1 and rangeTest == 1:
			simulation_report_f.write("Good simulation found %d!!!\n"%n)
			simulation_report_f.close()
			return n

		n+=1

	# if didn't find any good sim, return 0
	simulation_report_f.write("No Good Sim Found !!!\n")
	simulation_report_f.close()
	return 0

#----------------------------------------------------------------------------------------------------------#
def END_file(working_dir,string):
	f_end = open(working_dir+'/END.txt','a')
	f_end.write(string + '\n')
	return
#----------------------------------------------------------------------------------------------------------#
def getBestModel(model_data_dict):
	# This function receives the data read from ChromRes files for each model and returns the model
	#with the lowest AIC value.
	Best_Model_LowestAIC = 'NA'
	LowestAic = 99999999
	for model_key in model_data_dict.keys():
		if float(model_data_dict[model_key][2]) < LowestAic:
			Best_Model_LowestAIC = model_key
			LowestAic = float(model_data_dict[model_key][2])

	return Best_Model_LowestAIC,LowestAic

#----------------------------------------------------------------------------------------------------------#
def get_allFreq(line):
	#This function extract the chromCount frequenceis from the chromRes line of frequencies:
	#it returns a dictionary: keys are the chrom count and val is the freq val

	#F[21] = 9.12831e-06     F[22] = 9.13768e-06     F[23] = 9.16383e-06
	freq_dict=dict()
	line_split = line.split()
	for item in line_split:
		number = re.search("F\[(\d+)\]", item)
		freq = re.search("=(-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?)",item)
		if number:
			freq_dict[number.group(1)] = freq.group(1)
	return freq_dict
#----------------------------------------------------------------------------------------------------------#
def summarize_results(working_dir,tree_dir,models_list,param_dict): # check perl .pm parse_chromevol_results
	#For each tree directory we need to collect all chrom_res file results for analysis:

	#For debug:
	debug_f = open(tree_dir+'/DEBUG_file.txt','w')
	debug_f.write(working_dir)
	for model in models_list:
		debug_f.write(model)

	model_data_dict = dict()
	#tree_dir = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/SERVER_chromEvol/ChromEvol_Tree_1'
	#models_list = ['BASE_NUM','CONST_RATE']
	for model in models_list:
		chrom_res_f = tree_dir + '/' + model + '/chromEvol.res'
		#chrom_res_f = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/SERVER_chromEvol/ChromEvol_Tree_1/CONST_RATE/chromEvol.res'   #For debug'
		try:
			with open(chrom_res_f, 'r') as chrom_res_f:
				freqs = tree_dir + '/' + model + '/root_freq'
				gain = 'na';             gainL = 'na';            loss = 'na'
				lossL = 'na';            dupl = 'na';            demi = 'na'
				baseNum = 'na';          baseNumR = 'na'
				for line in chrom_res_f: # lnLik; AIC; totTrLen; rootFreq;
					line=line.strip()
					if 'LogLikelihood' in line:
						lnLik = line.split(' = ')[1]
						debug_f.write(lnLik+'\n')
					if 'AIC (Akaike information criterion)' in line:
						AIC = line.split(' = ')[1]
						debug_f.write(AIC + '\n')
					if 'GAIN_CONST' in line:
						gain = line.split()[1]
					if 'GAIN_LINEAR' in line:
						gainL = line.split()[1]
					if 'LOSS_CONST' in line:
						loss = line.split()[1]
					if 'LOSS_LINEAR' in line:
						lossL = line.split()[1]
					if 'DUPL' in line:
						dupl = line.split()[1]
					if 'HALF_DUPL' in line:
						demi = line.split()[1]
					if 'BASE_NUMBER_R' in line:
						baseNumR = line.split()[1]
					if 'BASE_NUMBER' in line:
						baseNum = line.split()[1]
						debug_f.write(baseNum + '\n')
					if '#total tree length' in line:
						totTrLen = line.split(' = ')[1]
					if '#max chromosome in data:' in line:
						maxChrNum = int(line.split(':')[1].strip())
					if 'F[' in line:
						create_freq_file(line, freqs, maxChrNum,param_dict,working_dir)
						rootFreq = get_allFreq(line)
						debug_f.write(str(rootFreq) + '\n')
				#Save dat for each model:
				model_data_dict[model] = [model, lnLik, AIC, totTrLen, gain, gainL, loss, lossL, dupl, demi, baseNum,
										  baseNumR, rootFreq]
				chrom_res_f.close()
		except IOError:
			return "Failed to open chrom_res file at: %s" %(tree_dir + '/' + model)

	#For debug:
	with open(tree_dir+'/DEBUG_file.txt','a') as debug_f:
		debug_f.write(working_dir)
		for model in models_list:
			debug_f.write(model_data_dict[model][0])
			debug_f.write(model_data_dict[model][1])
			debug_f.write(model_data_dict[model][2])
			debug_f.write(model_data_dict[model][3])

	BestModel,LowestAic = getBestModel(model_data_dict)
	BEST_MODEL_TREE_dict[tree_dir]=BestModel
	with open(tree_dir + '/ChosenModel_LowAIC.txt', 'w') as chosenModel_f:
		chosenModel_f.write(BestModel)
		if 'NO_DUPL' in BestModel:
			f_nodupl = open(tree_dir + '/NO_DUPL','w')
			f_nodupl.close()
	chosenModel_f.close()

	#Create result summary file:
	selected_models = ''
	dist_from_best_dict = dict()	#distance of AIC of each model from Best model
	with open(tree_dir+'/result_sum.csv','w') as sum_f:
		sum_f.write("model,LogLikelihood, AIC,AIC-LowestAIC,TreeLength,gain,gainL,loss,lossL,dupl,demi,baseNum,baseNumR\n")
		for model in model_data_dict.keys():
			if DEBUG_FLAG==1: print(model_data_dict[model])
			model_str=model_data_dict[model][0]+','+model_data_dict[model][1]+','+model_data_dict[model][2]+','+ \
					   str(float(model_data_dict[model][2])-LowestAic) + ','+model_data_dict[model][3]+','+ \
					   model_data_dict[model][4]+','+model_data_dict[model][5] + ','+model_data_dict[model][6]+','+ \
					   model_data_dict[model][7]+','+model_data_dict[model][8]+','+model_data_dict[model][9] + ','+ \
					   model_data_dict[model][10]+','+model_data_dict[model][11]
			sum_f.write(model_str+'\n')
			dist_from_best_dict[model] = format(float(model_data_dict[model][2]) - LowestAic,'.4f')
			if (float(model_data_dict[model][2])-LowestAic) == 0:
				selected_models+='%s,' %model
		#($tmp[0], $tmp[1], $tmp[2], $tmp[2] - $lowestAIC, $tmp[3], @tmp[4 .. 11])
	sum_f.close()
	shutil.copyfile(tree_dir+'/result_sum.csv', tree_dir+'/result_sum.txt')
	tree_number = re.search("ChromEvol_Tree_(\d+)$", tree_dir).group(1)
	#Check cases in which more than 1 model was selected (e.g CONST and CONST_DEMI)
	selected_models = selected_models[:-1] #remove final comma
	if len(selected_models.split(',')) > 1:
		selected_model = return_min_param_Model(selected_models)
	else:
		selected_model = selected_models.strip().replace(',','')


	with open(working_dir + '/Chromevol_selected_models.txt', 'a') as f_selectModels:
		f_selectModels.write("Selected model (AIC) for Tree# %s: %s\n" %(tree_number,selected_model))
	f_selectModels.close()
	return BestModel,dist_from_best_dict

#----------------------------------------------------------------------------------------------------------#
def send_job_cmd(cmd,job_dir,tree_model,queue_name):
	with open(job_dir + '/job_chrom.sh','w') as job_f:
		job_f.write("#!/bin/bash\n")
		job_f.write("#PBS -N chrom_%s\n" %tree_model)
		job_f.write("#PBS -r y\n")
		#job_f.write("#PBS -q lifesciweb\n")
		job_f.write("#PBS -q %s\n"%queue_name)
		job_f.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n")
		job_f.write("#PBS -e %s\n" %job_dir)
		job_f.write("#PBS -o %s\n" %job_dir)
		job_f.write("module load python/anaconda3-5.0.0\n")
		job_f.write("%s\n" %cmd)

	job_f.close()
	f_job_num = job_dir + '/Out_job_chrom.txt'
	with open(f_job_num,'w') as f_job:
		os.system("qsub %s" % job_dir + '/job_chrom.sh > ' + f_job_num)
		f_job.close()

	return read_job_id(f_job_num)
#----------------------------------------------------------------------------------------------------------#
def write_param_dict_to_file(param_dict,working_dir):

	dict_file = open(working_dir+'/param_dict.txt','w')
	for key in param_dict.keys():
		dict_file.write('%s,%s\n'%(key,param_dict[key]))
	return working_dir+'/param_dict.txt'


#----------------------------------------------------------------------------------------------------------#
def send_SIM_job_cmd(cmd,job_dir,Tree_NUM,queue_name):
	with open(job_dir + '/job_SIM_chrom.sh','w') as job_sim_f:
		job_sim_f.write("#!/bin/bash\n")
		job_sim_f.write("#PBS -N %d_chrom_SIM\n"%int(Tree_NUM))
		job_sim_f.write("#PBS -r y\n")
		job_sim_f.write("#PBS -q %s\n"%queue_name)
		job_sim_f.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n")
		job_sim_f.write("#PBS -e %s\n" %job_dir)
		job_sim_f.write("#PBS -o %s\n" %job_dir)
		job_sim_f.write("module load python/anaconda3-5.0.0\n")
		job_sim_f.write("%s\n" %cmd)

	job_sim_f.close()
	f_job_num = job_dir + '/Out_SIM_chrom.txt'
	with open(f_job_num,'w') as f_job:
		os.system("qsub %s" % job_dir + '/job_SIM_chrom.sh > ' + f_job_num)
		f_job.close()

	return read_job_id(f_job_num)

#----------------------------------------------------------------------------------------------------------#
def Check_if_ready(max_timeout,result_files):

	result_flag=0
	sleep_period = 2 * 60  # 2min
	poll_start_time = time.time()
	end_time = poll_start_time + (max_timeout * 3600)
	while time.time() < end_time:
		result_flag = 0
		print("Polling models ChromRes file ...time -> %s" % (time.time()))
		# Check if done:
		for res_file in result_files:
			if not os.path.exists(res_file):
				result_flag = 1
		if result_flag == 0:
			return 0
		else:
			time.sleep(sleep_period)
	return 1

#----------------------------------------------------------------------------------------------------------#
def Check_if_SIM_ready(max_timeout,working_tree_dirs,sim_summary):

	result_flag=0
	sleep_period = 2 * 60  # 2min
	poll_start_time = time.time()
	end_time = poll_start_time + (max_timeout * sleep_period)
	while time.time() < end_time:
		result_flag = 0
		# Check if done:
		for tree_dir in working_tree_dirs:
			if os.path.exists(tree_dir + '/NO_expectations'):
				sim_summary.write("%s sim completed - NO_expectations\n" % tree_dir)
				continue
			elif os.path.exists(tree_dir + '/NO_good_sim'):
				sim_summary.write("%s sim completed - NO_good_sim\n" % tree_dir)
				continue
			elif os.path.exists(tree_dir + '/NO_inference_step_complete'):
				sim_summary.write("%s sim completed - NO_inference_step_complete\n" % tree_dir)
				continue
			elif os.path.exists(tree_dir + '/inference_step_complete'):
				sim_summary.write("%s sim completed - inference_step_complete\n" % tree_dir)
				continue
			elif os.path.exists(tree_dir + '/NO_DUPL'):
				sim_summary.write("%s No duple - no simulation\n" % tree_dir)
				continue
			else:
				sim_summary.write("%s sim not completed - inference_step_complete\n" % tree_dir)
				result_flag = 1
		if result_flag == 0:
			return 0
		else:
			time.sleep(sleep_period)
	return 1

# ----------------------------------------------------------------------------------------------------------#
def Check_if_ready_Return(max_timeout, result_files):
	result_flag = 0
	sleep_period = 2 * 60  # 2min
	poll_start_time = time.time()
	end_time = poll_start_time + (max_timeout * 3600)
	while time.time() < end_time:
		result_flag = 0
		print("Polling models ChromRes file ...time -> %s" % (time.time()))
		# Check if done:
		for res_file in result_files:
			if not os.path.exists(res_file):
				result_flag = 1
				return
		if result_flag == 0:
			return 0
		else:
			time.sleep(sleep_period)
	return 1
#----------------------------------------------------------------------------------------------------------#
def check_and_add_param(line,param_name,param_dict):
	if param_name in line:
		param_val = line.strip().split(' ')[1]
		param_dict[param_name] = param_val
	return
#----------------------------------------------------------------------------------------------------------#
def get_models_from_file(params_dict):
	params_list = []
	if params_dict['_runModels']=="ALL":
		#params_list = ["CONST_RATE","CONST_RATE_DEMI","CONST_RATE_DEMI_EST","CONST_RATE_NO_DUPL",
		#			"LINEAR_RATE","LINEAR_RATE_DEMI","LINEAR_RATE_DEMI_EST","LINEAR_RATE_NO_DUPL",
		#			"BASE_NUM","BASE_NUM_DUPL"]
		params_list = ["CONST_RATE","CONST_RATE_DEMI","CONST_RATE_DEMI_EST","CONST_RATE_NO_DUPL",
						"BASE_NUM","BASE_NUM_DUPL"]
		return params_list
	elif params_dict['_runModels']=="CONST":
		params_list = ["CONST_RATE","CONST_RATE_DEMI","CONST_RATE_DEMI_EST","CONST_RATE_NO_DUPL",
					"BASE_NUM","BASE_NUM_DUPL"]
		return params_list
	elif params_dict['_runModels']=="LINEAR":
		params_list = ["LINEAR_RATE","LINEAR_RATE_DEMI","LINEAR_RATE_DEMI_EST","LINEAR_RATE_NO_DUPL"]
		return params_list
	else:
		split_line = params_dict['_runModels'].split(',')
		for item in split_line:
			if item != '':
				params_list.append(item)
		return params_list
#----------------------------------------------------------------------------------------------------------#
#def prepare_singleTree_run(ModelsList):
#
#    for model in ModelsList:
#
#
#
#    return
#
#----------------------------------------------------------------------------------------------------------#
def create_freq_file(line,freqs,emp_data,param_dict,working_dir):  # Written by Anna
# emp_data is a list of the min_chr_in_data,min_chr_allowed, max_... , max_...
	tmp = line
	# write the tmp to a file, with the aid of emp_data
	try:
		with open(freqs, "w+") as root_freq:
			text = tmp.split()
			first = re.search("F\[(\d+)\]", text[0])
			first = int(first.group(1))
			last = re.search("F\[(\d+)\]", text[-1])
			last = int(last.group(1))
			for i in range(1, first):
				root_freq.write("F[" + str(i) + "]=0\n")
			for i in range(0, last - first + 1):
				root_freq.write(text[i]+'\n')
			#for i in range(last + 1, emp_data[3] * 2 + 1):
			for i in range(last + 1, emp_data * 2 + 1):
				root_freq.write("F[" + str(i) + "]=0\n")
	except IOError:
		print("Can't open freqs file %s\n" % freqs)
		send_end_email(param_dict, working_dir,"FAIL")
		sys.exit()

	return
#----------------------------------------------------------------------------------------------------------#
def find_MinMax(working_dir,counts_file,lower_bound):
	minInData=0; maxInData=0; maxChromNum=0
	spc_list,counts_dict,counts_dict_all=get_spc_in_counts(working_dir,counts_file)
	key_list_toRemove=[]
	#remove counts lower than lower_bound:
	for key in counts_dict.keys():
		if counts_dict[key] == 'x':
			continue
		elif int(counts_dict[key]) < lower_bound:
			key_list_toRemove.append(key)
	for key_rem in key_list_toRemove:
		for key in counts_dict.keys():
			if key == key_rem:
				del counts_dict[key]
				break
	no_x_counts_list = {v: v for k, v in counts_dict.items() if v is not 'x'}
	minInData = int(min(no_x_counts_list))
	maxInData = int(max(no_x_counts_list))

	print_to_log(str(minInData), working_dir + '/Log.txt')
	print_to_log(str(maxInData), working_dir + '/Log.txt')

	maxChromNum = 2 * (int(maxInData) + 10)
	print("Line: 1091 - %s" % str(minInData))
	return minInData,maxInData,maxChromNum

#----------------------------------------------------------------------------------------------------------#
def create_dir_if_not_exists(dir,working_dir):
	if not os.path.exists(dir):
		os.makedirs(dir)
		print_to_log("Created Directory: %s" %dir, working_dir + '/Log.txt')

#----------------------------------------------------------------------------------------------------------#
def get_leaf_names(NewickTree):
	spc_names_on_Tree=[]
	tree = Phylo.read(NewickTree, "newick")
	for leaf in tree.get_terminals():
		spc_names_on_Tree.append(leaf.name)
	return spc_names_on_Tree

#----------------------------------------------------------------------------------------------------------#
def print_to_log(lineData, Log_f):

	with open (Log_f,'a') as f_log:
		f_log.write(lineData + '\n')
	f_log.close()
	return
#----------------------------------------------------------------------------------------------------------#
def init_param_dict():
	param_dict=dict()
	with open('/bioseq/chromEvol/power_PARAM_templates/DEFAULT_PARAM_VALUES','r') as def_f:
		for line in def_f:
			if line.strip():
				line_split = line.strip().split(':')
				param_dict[line_split[0]]=line_split[1]
	return param_dict

#----------------------------------------------------------------------------------------------------------#
def read_pip_control(working_dir,pip_file): #read_input

	# Set dictionary with default values:
	param_dict = init_param_dict()

	param_list = ['_dataFile', '_treesFile', '_topologyNum', '_simNum', '_runModels', '_plot', '_ploidyModelsOnly',
				  '_ploidyCallType', '_paramTemplates', '_outDir', '_name', '_logFile', '_fixRoot', '_excludeModels',
				  '_cpusNum', '_chromevolExe', '_baseNumber', '_baseNum']

	# Check for Model Adequacy parameters in params.txt file
	MA_ParamList = []
	rerun_flag = 'RegRen'
	with open(working_dir + '/params.txt', 'r') as params_f:
		for line in params_f:
			param_name = line.strip().split(':')[0]
			param_val = line.strip().split(':')[1]
			param_dict[param_name] = param_val
			if 'Origin' in line:
				rerun_flag = 'Rerun'
			if 'MA_' in line:
				MA_ParamList.append(param_name)

	if rerun_flag == 'Rerun':
		shutil.copyfile(param_dict['OriginJobDir'] + '/PIP_control', working_dir + '/PIP_control')
		shutil.copyfile(param_dict['OriginJobDir'] + '/TreesFile.txt', working_dir + '/TreesFile.txt')
		shutil.copyfile(param_dict['OriginJobDir'] + '/countsFile_edit', working_dir + '/countsFile_edit')
		shutil.copyfile(param_dict['OriginJobDir'] + '/selected_models.txt', working_dir + '/selected_models.txt')
		shutil.copyfile(param_dict['OriginJobDir'] + '/taxa_input_report.txt', working_dir + '/taxa_input_report.txt')

	with open(pip_file, 'r') as pip_f:
		for line in pip_f:
			for param_name in param_list:
				check_and_add_param(line, param_name, param_dict)

	# Edit counts file:
	tree_for_counts = open(working_dir + '/counts_Tree', 'w')
	with open(param_dict['_treesFile'], 'r') as trees_f:
		for line in trees_f:
			if line.strip():
				tree_for_counts.write(line.strip())
				break
		tree_for_counts.close()

	# edit_counts_file will also check for indiscrepancies between the counts and tree file the user entered.
	# In case of a mismatch we produce the corresponding not to the user to enable correction:
	JobID = os.path.basename(os.path.normpath(working_dir))
	flag_mismatch, spc_no_counts, spc_missing_onTree = edit_counts_file(working_dir,
																		working_dir + '/counts_Tree',
																		working_dir + '/countsFile',
																		working_dir + '/countsFile_edit', JobID)

	# if param_dict['_baseNumber'] == -999:  #User didn't give it as an input
	param_dict['_dataFile'] = param_dict['_dataFile'] + '_edit'
	Min_in_val, Max_in_val, MaxChrom = find_MinMax(working_dir, param_dict['_dataFile'], 4)
	if param_dict['_fixRoot'] == 0:
		param_dict['_maxChrNum'] = -1
	elif (Max_in_val > int(param_dict['_fixRoot'])):
		param_dict['_maxChrNum'] = -10
	else:
		param_dict['_maxChrNum'] = param_dict['_fixRoot'] + 1
	# ?????????????????????????????? Check with Lior/Itay
	# param_dict['_maxChrNum'] = -10  #1+(2*Max_in_val)
	param_dict['_minChrNum'] = -1
	param_dict['_baseNumber'] = Min_in_val
	if param_dict['_baseNumber'] != 0:
		param_dict['_bOptBaseNumber'] = 1
	param_dict['Min_in_val'] = Min_in_val
	param_dict['Max_in_val'] = Max_in_val

	for param in param_dict.keys():
		line_to_print = ('%s:%s' % (param, param_dict[param]))
		print_to_log(line_to_print, working_dir + '/Log.txt')

	return param_dict,MA_ParamList,rerun_flag,flag_mismatch, spc_no_counts, spc_missing_onTree
#----------------------------------------------------------------------------------------------------------#
def turn_list_to_DB_list(list_name):

	db_names_str=''
	for item in list_name:
		name_lowercase = str(item).lower()
		item_add="'"+name_lowercase+"',"
		db_names_str+=item_add
	final_str = db_names_str[:-1]
	return final_str

#----------------------------------------------------------------------------------------------------------#
def search_counts_in_CCDB(taxa_list,working_dir):

	spc_found_on_CCDB = []
	count_counts=0 #number of species for which there is data in CCDB
	conn_db = sqlite3.connect(CCDB_DATABASE)
	db_curser = conn_db.cursor()
	taxa_list_str = turn_list_to_DB_list(taxa_list)
	query_taxa_counts_db = "SELECT distinct resolved_name,median from resolved_chrom_counts where lower(resolved_name) in (%s)" % taxa_list_str
	f_report = open(working_dir+'/Species_Counts_match.txt','a')
	#f_report.write("Check for missing counts data in CCDB:")
	#f_report.write('%s\n' % query_taxa_counts_db)
	# query_genus_species_db = "SELECT id,genus,resolved_name_full,resolved_name,parsed_n,median,status from resolved_chrom_counts WHERE genus is '%s'" % genus
	# query_genus_species_db = "SELECT resolved_name FROM resolved_chrom_counts WHERE genus is '%s'" % genus
	db_curser.execute(query_taxa_counts_db)
	rows_genus_db = db_curser.fetchall()
	for line in rows_genus_db:
		specie_name = line[0]
		median_val = line[1]
		spc_found_on_CCDB.append(specie_name)
		f_report.write('>%s\n' % specie_name.replace(' ', '_'))
		f_report.write('%s\n' % median_val)
		count_counts += 1

	conn_db.close()
	return spc_found_on_CCDB
#----------------------------------------------------------------------------------------------------------#
def search_spc_in_NCBI(taxa_list,working_dir,rep_f):

	count_counts=0 #number of species for which there is data in NCBI
	spc_with_seqData=[]
	conn_db = sqlite3.connect(NCBI_DB_DATABASE)
	db_curser = conn_db.cursor()
	taxa_list_str = turn_list_to_DB_list(taxa_list)
	query_taxa_counts_db = "SELECT DISTINCT organism FROM Genbank_seqs_pln WHERE LOWER(organism) in (%s)" % taxa_list_str
	f_report = open(working_dir+'/Species_Counts_match.txt','a')
	f_report.write("Check for missing counts data in CCDB:")
	f_report.write('%s\n' % query_taxa_counts_db)
	# query_genus_species_db = "SELECT id,genus,resolved_name_full,resolved_name,parsed_n,median,status from resolved_chrom_counts WHERE genus is '%s'" % genus
	# query_genus_species_db = "SELECT resolved_name FROM resolved_chrom_counts WHERE genus is '%s'" % genus
	db_curser.execute(query_taxa_counts_db)
	rows_ncbi_db = db_curser.fetchall()
	for line in rows_ncbi_db:
		taxa_found = line[0]
		if taxa_found not in spc_with_seqData:
			spc_with_seqData.append(taxa_found)

	conn_db.close()
	return spc_with_seqData

#----------------------------------------------------------------------------------------------------------#
def edit_counts_file(working_dir,chosenTree, countsFile, outCounts, jobId):

	#Need to check 4 cases:
	#1. Species missing in counts file -> turn to X data
	#2. Species missing in the tree file -> remove from counts file
	#3. Species with missing counts that have data in CCDB
	#4. Species missing from the phylogeny with seq data at NCBI (usign OneTwoTree you can generate new phylogeny that will include those species)
	f_inputs_taxa_report = open(working_dir+'/taxa_input_report.txt','w')

	#Get species names from tree file
	spc_list_on_tree=get_leaf_names(chosenTree)
	#working_dir=''
	spc_list_in_counts,counts_dict,counts_all = get_spc_in_counts(working_dir,countsFile)
	flag_at_least_one=0
	species_without_counts = []
	species_with_counts = []
	#Edit counts file only if use asked to add 'x' for missing taxa:
	with open(outCounts,'w') as countsEdit_f:
		for spc_name in spc_list_on_tree:
			if spc_name in spc_list_in_counts:
				countsEdit_f.write('>%s\n' % spc_name)
				countsEdit_f.write('%s\n' % str(counts_dict[spc_name]))
				species_with_counts.append(spc_name)
			else:
				#countsEdit_f.write('>%s\n' % spc_name)
				#countsEdit_f.write('x\n')
				species_without_counts.append(spc_name.replace('_',' '))



	shutil.copyfile(countsFile,working_dir+'/counts_origin')
	shutil.copyfile(outCounts,countsFile)
	#check if species have counts on ccdb:
	if species_without_counts:
		spc_found_ccdb = search_counts_in_CCDB(species_without_counts, working_dir)
		f_inputs_taxa_report.write("The following species were found in your phylogeny but are missing in "
		"your counts file:<br>\n %s<br>\n" %(', '.join(map(str, species_without_counts))))
		f_inputs_taxa_report.write("<a href=\"chromEvol_results/" + jobId + "/RemoveSpc.txt\" target=\"_blank\">[Missing Spc],</a>")
		if spc_found_ccdb:
			f_inputs_taxa_report.write('Note that the following speices have counts data in '
			'<a href="http://ccdb.tau.ac.il/" target="_blank"> CCDB:</a><br>\n %s<br>\n' %
			(', '.join(map(str, spc_found_ccdb))))
		flag_at_least_one = 1
		with open(working_dir + '/Status.txt', 'w') as stats_f:
			stats_f.write("We found a mismatch between your phylogeny and counts file, see detailes below.")

	#Check if speices in counts file are missing from phylogeny:
	spc_missing_from_phylogeny=list(set(spc_list_in_counts)-set(spc_list_on_tree))
	if spc_missing_from_phylogeny:
		f_inputs_taxa_report.write("The following species were found in your counts file but are missing in "
		"your phylogeny file:<br>\n %s<br>\n" %(', '.join(map(str, spc_missing_from_phylogeny))))
		spc_missing_from_phylogeny_no_underscore = [x.replace('_',' ') for x in spc_missing_from_phylogeny]
		spc_found_ncbiDB = search_spc_in_NCBI(spc_missing_from_phylogeny_no_underscore, working_dir,f_inputs_taxa_report)
		if spc_found_ncbiDB:
			f_inputs_taxa_report.write('Note that the following speices have sequence data in '
			'<a href="https://www.ncbi.nlm.nih.gov/taxonomy" target="_blank"> NCBI:</a><br>\n %s. (You can create an '
			'updated phylogeny using <a href="http://onetwotree.tau.ac.il/" target="_blank"> OneTwoTree</a>).' %
			(', '.join(map(str, spc_found_ncbiDB))))
		# Add remark to correct the phylogeny tree to fit the counts file:
		flag_at_least_one = 1
		with open(working_dir + '/Status.txt', 'w') as stats_f:
			stats_f.write("We found a mismatch between your phylogeny and counts file, see detailes below.")

	if flag_at_least_one != 1:	#counts and phylogeny match:
		f_inputs_taxa_report.write("\n")

	f_inputs_taxa_report.close()
	return flag_at_least_one,species_without_counts,spc_missing_from_phylogeny
#----------------------------------------------------------------------------------------------------------#
def choose_random_trees(working_dir,treesFile, selected_TreeNum):

	trees_list =[]
	dirs_list = []
	dirs_tree_dict=dict()
	with open(treesFile,'r') as trees_f:
		for line in trees_f:
			if line.strip():
				trees_list.append(line.strip())
	total_trees_num = len(trees_list)
	if selected_TreeNum <= total_trees_num:
		print_to_log("%d Trees chosen out of %d" % (selected_TreeNum,total_trees_num),
					 working_dir + '/Log.txt')
		if total_trees_num == 1:
			chosen_trees = [0]
		else:
			chosen_trees = random.sample(range(1, total_trees_num), selected_TreeNum)
		index_num = 1
		for num in chosen_trees:
			Tree_dir = working_dir + '/ChromEvol_Tree_' + str(index_num)
			dirs_list.append(Tree_dir)
			Tree_file = Tree_dir + '/tree_' + str(index_num)
			dirs_tree_dict[Tree_dir]=Tree_file
			index_num+=1
			create_dir_if_not_exists(Tree_dir, working_dir)
			tree_f = open(Tree_file,'w')
			tree_f.write(trees_list[num])
			tree_f.close()
	else:
		print_string = "Not enough trees in input file: %d (chosen Trees %d)" % (total_trees_num,selected_TreeNum)
		END_file(working_dir, print_string)
		print_to_log("Not enough trees in input file: %d (chosen Trees %d)" % (total_trees_num,selected_TreeNum),
					 working_dir + '/Log.txt')

	return dirs_list,dirs_tree_dict
#----------------------------------------------------------------------------------------------------------#
def edit_param_file(working_dir,tree_dir,dirs_tree_dict,model_name,param_dict):
	#for model_item in params_list:
	param_local_dir = tree_dir #+ '/' + model_name
	param_local_file = param_local_dir + '/param_' + model_name + '.txt'
	if not os.path.exists(param_local_dir): os.makedirs(param_local_dir)
	shutil.copyfile(PARAMS_PATH + '/param_' + model_name, param_local_file)
	new_data = ''
	with open(param_local_file, 'r') as param_file:
		for line in param_file:
			new_line = line.replace('<OUT_DIR>', param_local_dir)
			new_line = new_line.replace('<CNT_FILE>', param_dict['_dataFile'])
			new_line = new_line.replace('<TREE_FILE>', dirs_tree_dict[tree_dir])
			new_line = new_line.replace('<MAX_CHR_NUM>', str(param_dict['_maxChrNum']))  # str(max_chrom_num))
			new_line = new_line.replace('<MIN_CHR_NUM>', str(param_dict['_minChrNum']))
			new_line = new_line.replace('<BASE_NUM>', str(param_dict['_baseNumber']))
			new_line = new_line.replace('<OPT_BASE_NUM>', str(param_dict['_bOptBaseNumber']))
			new_data += (new_line)
	with open(param_local_file, 'w') as param_file:
		param_file.write(new_data)

	return param_local_file

# ----------------------------------------------------------------------------------------------------------#
def get_spc_in_counts(working_dir,countsFile):

	spc_in_couns_list=[]
	counts_dict=dict()
	counts_dict_all=dict()
	flag_name = 0
	with open(countsFile,'r') as counts_f:
		for line in counts_f:
			print(line)
			if not line.strip():  continue   #To skip empty lines
			if '>' in line and flag_name == 0:
				flag_name = 1
				spc_name = line.strip().replace('>','')
				spc_in_couns_list.append(spc_name)
			else:
				if 'x' in line.strip() or 'X' in line.strip():
					counts_dict[spc_name] = line.strip().replace('X','x')
					counts_dict_all[spc_name] = line.strip().replace('X','x')
					none=1
				else:
					counts_dict[spc_name]=int(line.strip())
					counts_dict_all[spc_name]=int(line.strip())
				flag_name=0

		#In case this Genus has less than 5 counts:
		try:
			if len(spc_in_couns_list) < 5:
				open(working_dir + '/LessThanFiveWithCounts.txt', 'w')
				raise Exception
		except Exception:
			with open(working_dir + '/Status.txt', 'w') as stats_f:
				stats_f.write("This genus has less than 5 counts.")
			with open(working_dir + '/END.txt', 'a') as stats_f:
				stats_f.write("LessThanFiveWithCounts")
			return ValueError


	return spc_in_couns_list,counts_dict,counts_dict_all

# ----------------------------------------------------------------------------------------------------------#
def creat_models_and_params(working_dir,tree_dir,dirs_tree_dict,models_list,param_dict,queue_name):

	result_files=[]
	for model in models_list:
		model_dir = tree_dir + '/' + model
		print_to_log("Create model dir: %s" % model_dir, working_dir + '/Log.txt')
		create_dir_if_not_exists(model_dir,working_dir)
		param_file = edit_param_file(working_dir,tree_dir,dirs_tree_dict,model,param_dict)
		#Run Chromevol:
		cmd = '%s %s > %s/chrom_out.txt' % (param_dict['_chromevolExe'], param_file,model_dir)
		tree_model = model
		job_id = send_job_cmd(cmd, model_dir, tree_model,queue_name)
		print_to_log(job_id,working_dir + '/Log.txt')
		#print(cmd)
		#os.system(cmd)
		result_files.append(tree_dir + '/' + model + '/chromEvol.res')

	return result_files

def convertToNewModelNames(oldName): # josef
    if oldName=="CONST_RATE_NO_DUPL":
        return "Dys";
    if oldName=="CONST_RATE":
        return "DysDup";
    if oldName=="CONST_RATE_DEMI":
        return "DysDupDem=";
    if oldName=="CONST_RATE_DEMI_EST":
        return "DysDupDem";
    if oldName=="BASE_NUM":
        return "DysBnum";
    if oldName=="BASE_NUM_DUPL":
        return "DysBnumDup";
   
# ----------------------------------------------------------------------------------------------------------#
# ----------------------------------------------------------------------------------------------------------#
#
#
#													M A I N													#
#
#
# ----------------------------------------------------------------------------------------------------------#
# ----------------------------------------------------------------------------------------------------------#


#working_dir = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/SERVER_chromEvol/' #HEAD Directory / Job Dir
Log_dir = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/SERVER_chromEvol/' #HEAD Directory / Job Dir

#pip_control = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/SERVER_chromEvol/PIP_control'
#selected_models = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/SERVER_chromEvol/selected_models.txt'


working_dir = sys.argv[1]
pip_control = sys.argv[2]
queue_name = "lifesciweb"
selected_TreeNum = TREES_NUM #int(sys.argv[3]) #5  #number of trees to run on

try:
	#working_dir = working_dir + run_num
	create_dir_if_not_exists(working_dir,working_dir)

	#Check if chromevol results exists to skip to model adequacy stage:
	results_flag = 0
	for dirs in glob(working_dir):
		list_dirs = glob(working_dir + '/*/')
		for sub_dir in list_dirs:
			if 'ChromEvol_Tree' in sub_dir:
				#check if result file in place:
				if os.path.exists(sub_dir + 'result_sum.csv'):
					results_flag+=1

	results_lists = []
	param_dict,MA_param_list,Rerun_flag,flag_mismatch, spc_no_counts, spc_missing_onTree = read_pip_control(working_dir,pip_control)
	#Prepare email notifications cmds:
	if param_dict['jobTitle'] != 'daily test':
		send_start_email(param_dict,working_dir)

	Log_f = working_dir + '/Log.txt'
	#Set ploidy and model adeq params:
	if 'ploidy_ON' not in param_dict.keys():
		ploidy_flag = 'OFF'
	elif (param_dict['ploidy_ON'] == 'On'):
		ploidy_flag = 'ON'
	else:
		ploidy_flag = 'OFF'
	if (param_dict['AdeqTest'] == 'YES'):
		model_adq_flag = 'ON'
	else:
		model_adq_flag = 'OFF'
	print_to_log("Ploidy:%s,Adequacy:%s\n"%(ploidy_flag,model_adq_flag),working_dir + '/Log.txt')
	with open(working_dir+'/ma_ploidy_flags.txt','w') as ma_ploidy_f:
		ma_ploidy_f.write('%s,%s,'%(ploidy_flag,model_adq_flag))
	ma_ploidy_f.close()

	if (Rerun_flag == 'RegRun'):
		Flag1 = 'None'
	elif (Rerun_flag == 'Rerun'):
		with open(working_dir + '/Status.txt', 'w') as stats_f:
			stats_f.write("Rerun job %s"%param_dict['OriginJobID'])
		copy_dir_to_another_dir(param_dict['OriginJobDir'] +'/ChromEvol_Tree_1', working_dir, Log_f)


	# ----------------------------------------------------------------------------------------------------------#
	# ---------------------------------       Checking the input files         ---------------------------------#
	# ----------------------------------------------------------------------------------------------------------#
	if (Rerun_flag != 'Rerun'):
		with open(working_dir + '/Status.txt', 'w') as stats_f:
			stats_f.write("Running Chromevol")

	BestModels_perTree_dict = dict()
	bestModelsPerTree_f = open(working_dir+'/BestModel_per_Tree.txt','w')
	# ----------------------------------------------------------------------------------------------------------#
	if (Rerun_flag != 'Rerun'):
		#Edit counts file:
		tree_for_counts = open(working_dir + '/counts_Tree','w')
		with open(param_dict['_treesFile'], 'r') as trees_f:
			for line in trees_f:
				if line.strip():
					tree_for_counts.write(line.strip())
					break
			tree_for_counts.close()

		# edit_counts_file will also check for indiscrepancies between the counts and tree file the user entered.
		# In case of a mismatch we produce the corresponding not to the user to enable correction:
		JobID = os.path.basename(os.path.normpath(working_dir))
		flag_mismatch,spc_no_counts,spc_missing_onTree =  edit_counts_file(working_dir,working_dir + '/counts_Tree', working_dir + '/countsFile', working_dir + '/countsFile_edit',JobID)
		remove_spc_f = working_dir + '/RemoveSpc.txt'
		with open(remove_spc_f,'w') as rem_spc_f:
			for taxa in spc_no_counts:
				rem_spc_f.write(taxa.replace(' ','_') + '\n')

		#if (param_dict['PruneCounts'] == 'Xcounts'):
		if (param_dict['PruneCounts'] == 'PruneTree') and os.stat(remove_spc_f).st_size != 0:
			#add_x_counts_for_missing_taxa() #working_dir, tree_file, remove_species_list
			Rcmd_pruneTree = ("R CMD BATCH '--args working_dir=" + '"' + working_dir + '"' +
								 " trees_file=" + '"' + working_dir +'/TreesFile.txt' + '"' + " remove_species_list=" + '"' + remove_spc_f + '"' +
								 " outfile=" + '"' +working_dir +'/TreesFile_R.txt' + '"' + "' " + CHROMEVOL_CODE_PATH + "prune.trees.R")
			print_to_log(Rcmd_pruneTree,Log_f)
			stats_r = os.system(Rcmd_pruneTree)
			if stats_r != 0:
				print_to_log("Prune R status is %s" % stats_r, Log_f)
				send_end_email(param_dict,working_dir,"FAIL")
				sys.exit()
			else:
				shutil.copyfile(working_dir +'/TreesFile.txt',working_dir +'/TreesFile_Origin.txt')
				shutil.copyfile(working_dir +'/TreesFile_R.txt',working_dir +'/TreesFile.txt')
		#if flag_mismatch == 1:
		#	with open(working_dir + '/Status.txt', 'w') as stats_f:
		#		stats_f.write("Check file Species_Counts_match.txt")
		#	sys.exit()


		models_list = get_models_from_file(param_dict)
		working_tree_dirs, dirs_tree_dict = choose_random_trees(working_dir, working_dir + '/TreesFile.txt',
																selected_TreeNum)
		#param_dict['_dataFile'] = param_dict['_dataFile']
		print_to_log("Number of trees to run on: %d" %len(working_tree_dirs), working_dir + '/Log.txt')

		if os.path.exists(working_dir+'/Over_50prec_phylogenies.txt'):
			create_ploidy(working_tree_dirs,working_dir)
			sys.exit()

		for tree_dir in working_tree_dirs:
			print_to_log("Create chromevol tree path: %s" % tree_dir, working_dir + '/Log.txt')
			#Set the models for each tree:
			#models_list = get_models_from_file(param_dict)
			print_to_log("models selected are: %s" % models_list, working_dir + '/Log.txt')
			results_lists.append(creat_models_and_params(working_dir,tree_dir,dirs_tree_dict,models_list,param_dict,queue_name))

		counter_time = 0
		while counter_time < 1000:
			for list_res in results_lists:
				if Check_if_ready(4,list_res) == 0:
					print_to_log("All models for this list are done !!!", working_dir + '/Log.txt')
				else:
					print_to_log("Some models are still missing (counter time = %d)" %counter_time, working_dir + '/Log.txt')
			counter_time+=1

		for tree_dir in working_tree_dirs:
			for model_name in models_list:
				dir_path = tree_dir + '/' + model_name + '/'
				input_tree_ml = tree_dir + '/' + model_name + '/mlAncestors.tree'
				input_tree_post = tree_dir + '/' + model_name + '/posteriorAncestors.tree'
				counts_path = working_dir + '/countsFile_edit'
				convert_to_PhyD3(input_tree_ml, counts_path, 'None')
				convert_to_PhyD3(input_tree_post, counts_path, 'None')

		with open(working_dir + '/Status.txt', 'w') as stats_f:
			stats_f.write("Running Chromevol")

		#Prepare simulations:
		param_dict_file = write_param_dict_to_file(param_dict, working_dir)
		tree_BestModel_dict=dict()
		for tree_dir in working_tree_dirs:

			BestModel,AIC_per_model_dict = summarize_results(working_dir, tree_dir, models_list,param_dict)
			tree_BestModel_dict[tree_dir]=BestModel
			BestModels_perTree_dict[tree_dir]=BestModel
			bestModelsPerTree_f.write("%s,"%BestModel)
			treeFile = dirs_tree_dict[tree_dir]

			#------------------------------     SIMULATIONS     ------------------------------#
			if 'NO_DUPL' not in BestModel:	#results_flag != selected_TreeNum:
				cmd_simulation = 'python %s/simulation_test.py -m %s -td %s -pd %s -tree %s -wd %s'\
				%(CHROMEVOL_CODE_PATH,BestModel,tree_dir,param_dict_file,treeFile,working_dir)
				Tree_NUM = tree_dir.split('_')[-1].replace('/','')
				send_SIM_job_cmd(cmd_simulation, tree_dir, Tree_NUM, queue_name)

		sim_summary = open(working_dir+'/Simulation_summary.txt','w')
		counter_sim_time=0
		while counter_sim_time < 1000:
			if Check_if_SIM_ready(2,working_tree_dirs,sim_summary) == 0:
				print_to_log("All simulation completed", working_dir + '/Log.txt')
				break
			else:
				print_to_log("Some simulations are still running", working_dir + '/Log.txt')
			counter_sim_time+=1

		print_to_log("Model Adeq flag is %s"%model_adq_flag, working_dir + '/Log.txt')
		model_adq_list = []
		if model_adq_flag != 'OFF':
			for tree_dir in working_tree_dirs:
				model_adq_list = MA_Test(tree_dir, tree_BestModel_dict[tree_dir], model_adq_flag, param_dict)


	#------------------------------     Check if ReRun     ------------------------------#

	print_to_log("Rerun_flag is %s" % Rerun_flag, working_dir + '/Log.txt')
	if (Rerun_flag is 'Rerun'):
		model_adq_list = []
		tree_BestModel_dict=dict()
		for tree_dir in working_tree_dirs:
			BestModel, AIC_per_model_dict = summarize_results(working_dir, tree_dir, models_list)
			tree_BestModel_dict[tree_dir] = BestModel
			BestModels_perTree_dict[tree_dir] = BestModel
			bestModelsPerTree_f.write("%s," % BestModel)
		if (Rerun_flag == 'Rerun' and model_adq_flag == 'ON'):
			model_adq_list = []
			for tree_dir in working_tree_dirs:
				model_adq_list = MA_Test(tree_dir, tree_BestModel_dict[tree_dir], model_adq_flag, param_dict)

	# -------------- Create Results for WEB -----------------#
	# Final file:
	f_final = open(working_dir + '/END.txt', 'a')
	f_final.write("Chromevol run was completed\n")

	print_to_log("Creating results files:", working_dir + '/Log.txt')
	print_to_log("models_list :%s"%models_list, working_dir + '/Log.txt')
	print_to_log("model_adq_list :%s"%model_adq_list, working_dir + '/Log.txt')
	Create_results_for_web(working_dir,models_list,model_adq_list,AIC_per_model_dict)
	#Create_results_for_web(working_dir,model_adq_list)
	with open(working_dir+'/Status.txt','w') as stats_f:
		stats_f.write("Results files ready")

	# ------------- Convert trees for new Phy3 --------------#
	#for tree_dir in working_tree_dirs:
	#	for model_id in models_list:
	#		tree_f = tree_dir + '/' + model_id + '/mlAncestors.tree'
	#		tree_xml = tree_dir + model_id + '/mlAncestors_phylo'
	#		Phylo.convert(tree_f, 'newick', tree_xml, 'phyloxml')


	#-------------- Create Ploidy -----------------#
	NoDuple = 'No'	# if more than 50% selected models are NoDuple -> this flag will be set to Yes

	#summarize_models:
	completedPhylogenies = check_inferences(working_dir,working_tree_dirs,working_dir)
	print_to_log("len(completedPhylogenies) is %d\n" %len(completedPhylogenies), working_dir + '/Log.txt')
	if len(completedPhylogenies) == 0:
		with open(working_dir + '/FAIL.txt','w') as outcome_f:
			outcome_f.write('Completed Phylogenies returned 0 -> Error while running\n')
		outcome_f.close()
		send_end_email(param_dict, working_dir,"FAIL")
		sys.exit()
	else:
		if len(completedPhylogenies) <= selected_TreeNum:
			modelsSumFile = working_dir + '/models_summary'
			models_sum = summarize_models(working_tree_dirs, modelsSumFile)
			if (models_sum == 0):
				NoDupleMoreThan50 = open(working_dir + '/NoDuple_MoreThan50.txt','w')
				NoDupleMoreThan50.write('More than 50% NoDuple models\n')
				NoDuple = 'Yes'
				NoDupleMoreThan50.close()
			elif len(completedPhylogenies) >= selected_TreeNum/2:
				Over_50p_Trees = open(working_dir + '/Over_50prec_phylogenies.txt', 'w')
				Over_50p_Trees.write('completedTreesNums: %d\n' % len(completedPhylogenies))
				Over_50p_Trees.close()
			else:
				NotEnoughTrees = open(working_dir + '/Less_than_50_phylogenies.txt','w')
				NotEnoughTrees.write('completedTreesNums: %d\n' %len(completedPhylogenies))
				NotEnoughTrees.close()
		else:
			print_to_log("completedPhylogenies equals %d\n" % selected_TreeNum, working_dir + '/Log.txt')



	print_to_log("Ploidy Flag - %s, NoDuple is %s\n"%(ploidy_flag,NoDuple),working_dir+'/Log.txt')
	if ploidy_flag == 'ON' and NoDuple == 'No':
		print_to_log("Call create_ploidy with %s\n"%working_tree_dirs, working_dir + '/Log.txt')
		create_ploidy(working_tree_dirs,working_dir)

	#--------------------------------------------------------------------------------------------------------------------#
	#		converting trees to phyD3
	#--------------------------------------------------------------------------------------------------------------------#

	for tree_dir in working_tree_dirs:
		for model_name in models_list:
			#dir_path = tree_dir + '/' + BestModels_perTree_dict[tree_dir] + '/'
			dir_path = tree_dir + '/' + model_name + '/'
			#input_tree = tree_dir + '/' + BestModels_perTree_dict[tree_dir] + '/mlAncestors.tree'
			input_tree = tree_dir + '/' + model_name + '/mlAncestors.tree'
			counts_path = working_dir + '/countsFile_edit'
			if ploidy_flag == 'ON' and NoDuple == 'No':
				ploidy_path = working_dir + '/ploidy.csv'
				convert_to_PhyD3(input_tree, counts_path, ploidy_path)
			else:
				convert_to_PhyD3(input_tree, counts_path, 'None')

			# Create zip files of results for each model for sownload:
			#with zipfile.ZipFile(dir_path + '/%s_results.zip'%model_name, 'a') as myzip:
			#	myzip.write(dir_path + '/Results_Events.csv','Results_Events.csv')
			#	myzip.write(dir_path + '/Results_ExpNum.csv','Results_ExpNum.csv')
			#	myzip.write(dir_path + '/mlAncestors.tree','mlAncestors.tree')
			#	myzip.write(dir_path + '/posteriorAncestors.tree','posteriorAncestors.tree')
			#	myzip.write(dir_path + '/chromEvol.res','chromEvol.res')
			#myzip.close()

	with open(working_dir+'/Status.txt','w') as stats_f:
		stats_f.write("Completed")
	with open(working_dir+'/JOB_PASS.txt','w') as pass_f:
		pass_f.write('\n')

	#End of run
	send_end_email(param_dict,working_dir,"PASS")

except:
	with open(working_dir+'/Status.txt','a') as stats_f:
		stats_f.write("Error while running, please contact us for more details.")
	with open(working_dir+'/END.txt','a') as stats_f:
		stats_f.write("Except error line 1868")
	send_end_email(param_dict, working_dir, "FAIL")


