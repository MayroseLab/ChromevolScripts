import argparse, os
from ChromEvol_defs import *
import json

PARAMS_PATH = '/bioseq/chromEvol/power_PARAM_templates'
DEBUG_FLAG=0
#-------------------------------------------------------------------------------------------------------------------#
def print_to_log(lineData, Log_f):

	with open (Log_f,'a') as f_log:
		f_log.write(lineData + '\n')
	f_log.close()
	return
#-------------------------------------------------------------------------------------------------------------------#
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


#-------------------------------------------------------------------------------------------------------------------#
#This fuction will create 200 simulation dirs from which 1 will be selected:
def prepare_simulations(working_dir,tree_dir,freqFile, paramTemplate, nbrSims,maxChrNum, simDir, param_dict):

	#resSum = tree_dir + '/result_sum.csv'
	bestModel = read_val_from_file(tree_dir + '/ChosenModel_LowAIC.txt',1)
	treeNum = re.search("(\d+)$", tree_dir).group(1)
	treeFile = tree_dir + "/tree_" + treeNum
	count = 1
	baseNum = None
	minChrNum = 1       # MD q: ????? why 1 and not the min entered in counts?????
	simParamFilesList = []

	if (bestModel is "CONST_RATE_NO_DUPL" or bestModel is "LINEAR_RATE_NO_DUPL"):
		f_flag = open(tree_dir+ '/NO_DUPL','w')
		f_flag.close()
		no_duple_print = "NO_DUPL model chosen, no ploidy inference will be performed for this phylogeny\n"
		print_to_log(no_duple_print, working_dir + '/Log.txt')

	#Create 200 simulations (assuming nbrSims == 100)
	while count <= 2*nbrSims:
		simNum = count
		# get param values
		paramValues = ""
		resFile = tree_dir + '/' + bestModel + '/chromEvol.res'
		with open(resFile, 'r') as res_f:
			for line in res_f:
				line=line.strip()
				if 'LOSS_CONST' in line:
					paramValues += '_lossConstR ' + line.split('\t')[1] + '\n'
				if 'GAIN_CONST' in line:
					paramValues += '_gainConstR ' + line.split('\t')[1] + '\n'
				if 'DUPL' in line:
					paramValues += '_duplConstR ' + line.split('\t')[1] + '\n'
				if 'BASE_NUMBER_R' in line:
					paramValues += '_baseNumberR ' + line.split('\t')[1] + '\n'
				if 'BASE_NUMBER' in line and not 'BASE_NUMBER_R' in line:
					paramValues += '_baseNumber ' + line.split('\t')[1] + '\n'
					baseNum = line.split('\t')[1]
				if 'HALF_DUPL' in line:
					paramValues += '_demiPloidyR ' + line.split('\t')[1] + '\n'
				if 'LOSS_LINEAR' in line:
					paramValues += '_lossLinearR ' + line.split('\t')[1] + '\n'
				if 'GAIN_LINEAR' in line:
					paramValues += '_gainLinearR ' + line.split('\t')[1] + '\n'
				if 'total tree length =' in line:
					totalTreeLength = line.split('=')[1]

		if(baseNum):
			expFile = tree_dir + '/' + bestModel + '/expectations.txt'
			mlAncTree = tree_dir + '/' + bestModel + '/mlAncestors.tree'
			bntpv = create_bntpv(expFile, mlAncTree, int(baseNum),bestModel)
			if bntpv:
				print(bntpv)
				paramValues += "_baseTransitionProbs " + bntpv + '\n'
			else:
				baseNum = 1
				print("1 -> no bntpv")
				send_end_email(param_dict,working_dir,"FAIL")
				sys.exit()

		### parse parameter file for ChromEvol simulation
		#create 1 out of 2*nbrSims simulation directories:
		outSimDir = simDir + '/dir_sim_' + str(simNum)
		create_dir_if_not_exists(outSimDir)
		### print param to file
		outParamFile = simDir + '/param_sim_' + str(simNum)

		new_data = ''
		with open(paramTemplate, 'r') as param_sim_f:
			for line in param_sim_f:
				line=line.strip()
				if '<FREQ_FILE>' in line: new_line = line.replace('<FREQ_FILE>', freqFile) + '\n'
				elif '<MIN_CHR_NUM>' in line: new_line = line.replace('<MIN_CHR_NUM>', str(minChrNum)) + '\n'
				elif '<MAX_CHR_NUM>' in line: new_line = line.replace('<MAX_CHR_NUM>', str(maxChrNum)) + '\n'
				elif '<PARAMETERS>' in line: new_line = line.replace('<PARAMETERS>', paramValues)  # str(max_chrom_num))
				elif '<TREE_LENGTH>' in line: new_line = line.replace('<TREE_LENGTH>', totalTreeLength.strip()) + '\n'
				elif '<MAX_CHR_NUM_SIM>' in line: new_line = line.replace('<MAX_CHR_NUM_SIM>', str(maxChrNum)) + '\n'
				elif '<OUT_DIR>' in line: new_line = line.replace('<OUT_DIR>', outSimDir) + '\n'
				elif '<TREE_FILE>' in line: new_line = line.replace('<TREE_FILE>', treeFile) + '\n'
				else: new_line=line + '\n'
				print(new_line)
				new_data += (new_line)

		with open(outParamFile, 'w') as param_file:
			param_file.write(new_data)
		simParamFilesList.append(outParamFile)

		count+=1

	return simParamFilesList

#----------------------------------------------------------------------------------------------------------#
def create_dir_if_not_exists(dir):
	if not os.path.exists(dir):
		os.makedirs(dir)

#-------------------------------------------------------------------------------------------------------------------#
def read_param_dict_to_dict(param_dict_file):
	dict_file = open(param_dict_file, 'r')
	param_dict=dict()
	for line in dict_file:
		key = line.strip().split(',')[0]
		#if '_runModels' in key:
		#	data = line.strip().split(',')[1:]
		if 'definedModels' in key:
			data = line.strip().split(',', 1)[1]
		else:
			data = line.strip().split(',')[1]
		param_dict[key]=data

	return param_dict

#-------------------------------------------------------------------------------------------------------------------#
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
def parse_sim_results(realCountsFile,simCountsFile,outFile,working_dir):

	### Run on chosen simulations
	print(realCountsFile)
	spc_in_couns_list,refData_partial,refData = get_spc_in_counts(working_dir,realCountsFile)
	### add 'x'
	spc_in_couns_list,simData_partial,simData = get_spc_in_counts(working_dir,simCountsFile)
	# replace simulated count with 'x' if present as 'x' in original chromosome number data
	for key in simData.keys():
		if refData[key] == 'x':
			simData[key] = 'x'

	with open(outFile,'w') as out_f:
		for key in simData.keys():
			out_f.write('>%s\n' %key)
			out_f.write('%s\n' %str(simData[key]))
	out_f.close()
	return
#----------------------------------------------------------------------------------------------------------#
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
				if 'x' in line.strip():
					counts_dict[spc_name] = line.strip()
					counts_dict_all[spc_name] = line.strip()
					none=1
				else:
					counts_dict[spc_name]=int(line.strip())
					counts_dict_all[spc_name]=int(line.strip())
				flag_name=0

	return spc_in_couns_list,counts_dict,counts_dict_all
#----------------------------------------------------------------------------------------------------------#
def get_models_from_file(params_dict):
	params_list = []
	print(params_dict['definedModels'])
	data = json.loads(params_dict['definedModels'])
	counter = 0
	for d in data: 
		params_list.append(d['name'])
		counter += 1
	return params_list, data
	
#----------------------------------------------------------------------------------------------------------#
def prepare_initial_run( working_dir, params_dict, treeFile, countFile, outDir):
#def prepare_initial_run(working_dir,params_dict,treeFile,countFile,paramDir,runModels,excludeModels,baseNum,root,outDir): # For simulations
#working_dir,param_dict,treeFile, simDataFile, paramTemplateDir, bestModel, 0, 0, 0, simInferDir,dirs_tree_dict
	# sub prepare initial chromevol run
	# prepares the required files for the initial chromevol run

	modelsList, data = get_models_from_file(params_dict)
	paramFiles = []
	for model in modelsList:    ### parameter files
		print(model)
		print("-----------------------------------------------> prepare_initial_run")
		paramFiles.append('param_' + model)
		param_file = 'param_' + model
		### ChromEvol run parameters
		### NB: NOT EXHAUSTIVE LIST OF PARAMETERS
		### 	SEE <http://www.tau.ac.il/~itaymay/cp/chromEvol/> FOR DETAILS
		#maxChrNum = -1

		#minInData, maxInData, maxChromNum = find_MinMax(working_dir, param_dict['_dataFile'], 1)
		#maxChrNum = -10
		#minChrNum = -1

		model_i = getIndexByModel(data, model);
		edit_param_file_sim1(outDir, treeFile, model, params_dict, working_dir, data[model_i])
		#edit_param_file_sim(working_dir, outDir, treeFile, model, param_dict)
		#parseParamFile(param_f, outDir, treeFile, countFile)


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
def prep_params_for_sim(working_dir,counts_sim,param_dict): #read_input
	#This function will calculate the params for simulations:

	# Set dictionary with default values:
	param_sim_dict = dict()

	# ['_maxChrNum','_minChrNum','_baseNumber','_bOptBaseNumber']

	# if param_dict['_baseNumber'] == -999:  #User didn't give it as an input
	Min_in_val, Max_in_val, MaxChrom = find_MinMax(working_dir, counts_sim, 4)
	if param_dict['_fixRoot'] == 0:
		param_sim_dict['_maxChrNum'] = -1
	elif (Max_in_val > int(param_dict['_fixRoot'])):
		param_sim_dict['_maxChrNum'] = -10
	else:
		param_sim_dict['_maxChrNum'] = param_dict['_fixRoot'] + 1
	# ?????????????????????????????? Check with Lior/Itay
	# param_dict['_maxChrNum'] = -10  #1+(2*Max_in_val)
	param_sim_dict['_minChrNum'] = -1
	param_sim_dict['_baseNumber'] = Min_in_val
	if param_sim_dict['_baseNumber'] != 0:
		param_sim_dict['_bOptBaseNumber'] = 1
	param_sim_dict['Min_in_val'] = Min_in_val
	param_sim_dict['Max_in_val'] = Max_in_val
	for param in param_sim_dict.keys():
		line_to_print = ('SIM:%s:%s' % (param, param_sim_dict[param]))
		print_to_log(line_to_print, working_dir + '/Log.txt')

	return param_sim_dict

#----------------------------------------------------------------------------------------------------------#
def edit_param_file_sim(out_dir,tree_file,model_name,param_dict,working_dir):

	#for model_item in params_list:
	param_local_dir = out_dir #+ '/' + model_name
	param_local_file = param_local_dir + '/param_' + model_name + '.txt'
	counts_file = param_local_dir + '/sim_data.counts'
	if not os.path.exists(param_local_dir): os.makedirs(param_local_dir)
	shutil.copyfile(PARAMS_PATH + '/param_' + model_name, param_local_file)
	new_data = ''

	param_dict_sim=prep_params_for_sim(working_dir, counts_file, param_dict)

	with open(param_local_file, 'r') as param_file:
		for line in param_file:
			if '_simulationsNum' in line:
				new_line = line.replace('_simulationsNum 0', '_simulationsNum 100')
			new_line = line.replace('<OUT_DIR>', param_local_dir)
			new_line = new_line.replace('<CNT_FILE>', counts_file)
			new_line = new_line.replace('<TREE_FILE>', tree_file)
			new_line = new_line.replace('<MAX_CHR_NUM>', str(param_dict_sim['_maxChrNum']))  # str(max_chrom_num))
			new_line = new_line.replace('<MIN_CHR_NUM>', str(param_dict_sim['_minChrNum']))
			new_line = new_line.replace('<BASE_NUM>', str(param_dict_sim['_baseNumber']))
			new_line = new_line.replace('<OPT_BASE_NUM>', str(param_dict_sim['_bOptBaseNumber']))
			new_data += (new_line)
	with open(param_local_file, 'w') as param_file:
		param_file.write(new_data)

	return param_local_file

#def edit_param_file1(working_dir,tree_dir,dirs_tree_dict,model_name,model_dict,param_dict):
def edit_param_file_sim1(out_dir, tree_file, model_name, param_dict, working_dir, model_dict):
	#for model_item in params_list:
	param_local_dir = out_dir #+ '/' + model_name
	param_local_file = param_local_dir + '/param_' + model_name + '.txt'
	if not os.path.exists(param_local_dir): os.makedirs(param_local_dir)
	shutil.copyfile(PARAMS_PATH + '/param_NEW_MODELS', param_local_file)
	counts_file = param_local_dir + '/sim_data.counts'
	new_data = ''
	
	print_to_log("calling prep_params_for_sim", working_dir + '/Log.txt')
	param_dict_sim=prep_params_for_sim(working_dir, counts_file, param_dict)
	
	with open(param_local_file, 'r') as param_file:
		for line in param_file:
			if '_simulationsNum' in line:
				new_line = line.replace('_simulationsNum 0', '_simulationsNum 100')
			new_line = line.replace('<OUT_DIR>', f'{param_local_dir}/{model_name}/')
			new_line = new_line.replace('<CNT_FILE>', counts_file)
			new_line = new_line.replace('<TREE_FILE>', tree_file)
			new_line = new_line.replace('<MAX_CHR_NUM>', str(param_dict_sim['_maxChrNum']))  # str(max_chrom_num))
			new_line = new_line.replace('<MIN_CHR_NUM>', str(param_dict_sim['_minChrNum']))
			if model_dict['_baseNumberR'] == 1: 
				logVal = 10
			else: 
				logVal = 6
			new_line = new_line.replace('<LOG_VAL>', str(logVal))
			new_data += (new_line)
			
	with open(param_local_file, 'w') as param_file:
		param_file.write(new_data)
		for key, value in model_dict.items():
			if value != 0 and key != "name": 
				param_file.write(f'{key} {value}\n')

		# add lines for _baseNumber, _bOptBaseNumber, _minBaseTransition according to model
		if model_dict['_baseNumberR'] == 1:
			param_file.write(f'_baseNumber {param_dict["_baseNumber"]}\n')
			if param_dict['_baseNumber'] != 0:
				param_file.write(f'_bOptBaseNumber {param_dict["_bOptBaseNumber"]}\n')
			param_file.write('_minBaseTransition 0\n')
			
	param_file.close()
	
	return param_local_file
	
#-----------------------------------------     Simulation for Chromevol    --------------------------------------#

# process input from command line
parser = argparse.ArgumentParser(
	description='Extract counts data from CCDB and Ploidy inference from Ploidb reults')
parser.add_argument('--best_model', '-m', help='Name of best model',required=True)
parser.add_argument('--tree_dir', '-td', help='directory of chromevol tree', required=True)
parser.add_argument('--param_f', '-pd', help='parameters dictionary file', required=True)
parser.add_argument('--tree_file', '-tree', help='tree file', required=True)
parser.add_argument('--work_dir', '-wd', help='working directory', required=True)

args = parser.parse_args()
BestModel = args.best_model
tree_dir = args.tree_dir
param_dict_file = args.param_f
working_dir = args.work_dir
treeFile = args.tree_file


param_dict=read_param_dict_to_dict(param_dict_file)
#if 'NO_DUPL' not in BestModel:  # results_flag != selected_TreeNum:
print("Summary is ready at %s !!!" % tree_dir)

# tree_dir,freqFile, paramTemplate, nbrSims,maxChrNum,chromNumRange, simDir, fake
freqFile = tree_dir + '/' + BestModel + '/root_freq'
paramTemplate = PARAMS_PATH + '/param_SIM'
nbrSims = 50
maxChrNum = (int(param_dict['Max_in_val']) + 10) * 2
chromNumRange = int(param_dict['Min_in_val']) - int(param_dict['Min_in_val'])
simDir = tree_dir + '/simulation'
create_dir_if_not_exists(simDir)
fake = 0
simParamFilesList = prepare_simulations(working_dir, tree_dir, freqFile, paramTemplate, nbrSims, maxChrNum, simDir,
										param_dict)
print(simParamFilesList)
for sim_file in simParamFilesList:
	sim_cmd = param_dict['_chromevolExe'] + ' ' + sim_file
	os.system(sim_cmd)
goodSim = test_simulations(nbrSims, simDir)  # , maxInData, minInData, simulationDir)
if goodSim == 0:
	f_no_sim = open(tree_dir+'/NO_good_sim','w')
	f_no_sim.close()

# if a good simulation was found
if (goodSim):
	## remove all but one good simulation
	for n in range(1, nbrSims * 2 + 1, 1):
		print(n)
		print(goodSim)
		print('%s/dir_sim_%s' % (simDir, str(n)))
		print('%s/param_sim_%s' % (simDir, str(n)))
		if n != goodSim:
			os.system('rm -rf %s/dir_sim_%s' % (simDir, str(n)))
			os.system('rm -rf %s/param_sim_%s' % (simDir, str(n)))

if str(goodSim) != '0':
	## prepare infer from sim
	goodSimDir = '%s/dir_sim_%s/0' % (simDir, str(goodSim))
	simCountsFile = goodSimDir + '/simCounts.txt'
	simInferDir = simDir + '/sim_infer'
	create_dir_if_not_exists(simInferDir)
	simDataFile = simInferDir + '/sim_data.counts'  # same as simCounts, but with missing species replaced with x
	parse_sim_results(param_dict['_dataFile'], simCountsFile, simDataFile, working_dir)
	#treeFile = dirs_tree_dict[tree_dir]
	#working_dir,params_dict,treeFile,outDir
	prepare_initial_run(working_dir, param_dict, treeFile, simDataFile, simInferDir)
	#prepare_initial_run(working_dir, param_dict, treeFile, simDataFile, PARAMS_PATH, BestModel, 0, 0, 0, simInferDir)
	#prepare_initial_run(working_dir, param_dict, treeFile, simInferDir)

	## run inference from sim
	runFile = simInferDir + '/param_' + BestModel + '.txt'
	print(param_dict['_chromevolExe'])
	print(runFile)
	sim_cmd = param_dict['_chromevolExe'] + ' ' + runFile
	os.system(sim_cmd)
	# Need to add polling function that will enable other jobs to continue while this one is still running:
	f_infer_step = open(tree_dir + '/inference_step_complete', 'w')
	f_infer_step.close()
