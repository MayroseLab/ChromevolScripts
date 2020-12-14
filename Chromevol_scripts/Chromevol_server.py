import sys
import shutil
import os
import time

CHROMEVOL_EXE = '/bioseq/chromEvol/chromEvol.exe'
PARAMS_PATH = '/bioseq/chromEvol/power_PARAM_templates'
PARAM_NAMES = ['SIM','SIM_NoDUPL','BASE_NUM_DUPL','BASE_NUM','LINEAR_RATE_NO_DUPL','LINEAR_RATE_DEMI_EST',
				'LINEAR_RATE_DEMI','LINEAR_RATE','CONST_RATE_NO_DUPL','CONST_RATE_DEMI_EST','CONST_RATE_DEMI',
				'CONST_RATE']


###########################################################################################
# Polling Jmodel jobs to check if they are done:
def polling_for_resultFile(max_timeout,result_file):

	#max_timeout = 1 # In hours
	sleep_period = 10*60    #10min
	poll_start_time = time.time()
	end_time = poll_start_time + (max_timeout * 3600)
	while time.time() < end_time:
		print("Polling %s ...time -> %s" % (result_file,time.time()))
		#Check if done:
		if os.path.exists(result_file):
			return 0
		else:
			time.sleep(sleep_period)
	return 1
###########################################################################################
# Polling Jmodel jobs to check if they are done:
def create_chrom_results_for_site(path,chrom_res_file):

	f_chrom_res = open(chrom_res_file,'r')
	with open(path + "/web_file_res.txt", 'w') as web_result:
		for line in f_chrom_res:
			if "#Input" in line or "F[" in line:
				continue
			else:
				web_result.write(line)
	web_result.close()
	return
###########################################################################################

#This file will create the files/paths/params required for a Chromevol run
#
# ../chromEvol.exe ../param_MODEL_NAME.txt, example:
# /bioseq/chromEvol/chromEvol.exe .../param_CONST_RATE_DEMI.txt

job_dir = sys.argv[1]
tree_file = sys.argv[2]
counts_file = sys.argv[3]

#check which models were selected by the user:
parmas_file = job_dir + '/selected_models.txt'
params_list = []
with open (parmas_file,'r') as params_f:
	for line in params_f:
		split_line = line.strip().split(',')
		for item in split_line:
			params_list.append(item)

print(params_list)

#Need to calc min max chrom_num:
max_chrom_num = -1
min_chrom_num = 999999
with open (counts_file, 'r') as counts_f:
	for line in counts_f:
		print(line)
		line = line.strip()
		if line != '' and '>' not in line:
			if int(line) > max_chrom_num:
				max_chrom_num = int(line)
			if int(line) < min_chrom_num:
				min_chrom_num = int(line)

max_chrom_num = -1 * max_chrom_num      # WHY should be -10 with relation to root !!!!!!!!!!!!!!!!!!!!
print(max_chrom_num)
print(min_chrom_num)

if max_chrom_num == -1:
	print("max chrom is -1 !!!")
	sys.exit()
#min_chrom_num = -1

#Copy param template to directory:
for model_item in params_list:
	param_local_dir = job_dir + '/' + model_item
	param_local_file = param_local_dir + '/param_' + model_item + '.txt'
	if not os.path.exists(param_local_dir): os.makedirs(param_local_dir)
	shutil.copyfile(PARAMS_PATH + '/param_' + model_item,param_local_file)
	new_data=''
	with open (param_local_file, 'r') as param_file:
		for line in param_file:
			new_line = line.replace('<OUT_DIR>',job_dir)
			new_line = new_line.replace('<CNT_FILE>',counts_file)
			new_line = new_line.replace('<TREE_FILE>',tree_file)
			new_line = new_line.replace('<MAX_CHR_NUM>','-10') #str(max_chrom_num))
			new_line = new_line.replace('<MIN_CHR_NUM>',str(min_chrom_num))
			new_line = new_line.replace('<BASE_NUM>',str(min_chrom_num))
			new_line = new_line.replace('<OPT_BASE_NUM>','1')
			new_data+=(new_line)
	with open (param_local_file, 'w') as param_file:
		param_file.write(new_data)

	# with open param_file
	cmd = '%s %s' %(CHROMEVOL_EXE,param_local_file)
	os.system(cmd)

#for model_item in params_list:
	model_path = job_dir + '/' + model_item + '/'
	with open(model_path + 'FINAL_RESULTS_FILE.txt','w') as final_f:
		final_f.write("asfjhalsjfhajfh")
	final_f.close()
	result_file = model_path + 'chromEvol.res'
	create_chrom_results_for_site(model_path, result_file)
	if polling_for_resultFile(1,result_file) == 0:
		shutil.copyfile(model_path + 'posteriorAncestors.tree',model_path + 'posteriorAncestors')
		shutil.copyfile(model_path + 'mlAncestors.tree',model_path + 'mlAncestors')
		shutil.copyfile(model_path + 'exp.tree',model_path + 'exp')
		shutil.copyfile(model_path + 'allNodes.tree',model_path + 'allNodes')


