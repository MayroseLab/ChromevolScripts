import random
import os
import sys
import shutil
import time
import re
from ete3 import Tree
from Bio import Phylo
import numpy as np



#----------------------------------------------------------------------------------------------------------#
def largest_value_hash(hash):
	# sub largest value in hash
	# Receives a hash (ref) with numeric values, and returns the key with the highest value
	big = 0
	for key in hash.keys():
		if hash[key] > big:
			big = hash[key]
			ret_key = key
	return ret_key
#----------------------------------------------------------------------------------------------------------#
def majority_in_arr(arr):
# sub majority in array
# receives an array (ref) and returns the most frequent element and it's fraction of the array
	majHash = dict()
	for n in arr:
		if n in majHash.keys():
			majHash[n]+=1
		else:
			majHash[n]=1

	mostFreqKey = largest_value_hash(majHash)
	freq = majHash[mostFreqKey]
	majority = freq/len(arr)

	return mostFreqKey,majority
#----------------------------------------------------------------------------------------------------------#
def MCC(TP, TN, FP, FN):
	# sub MCC
	# calculates Mathew's Correlation Coefficient
	MCC = 0
	# if MCC denominator is zero, set MCC value to zero
	if ((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) != 0):
		MCC = ((TP*TN) - (FP*FN)) / ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5

	return MCC
#----------------------------------------------------------------------------------------------------------#
def read_val_from_file(file_path,num_param_per_line):
	with open(file_path,'r') as r_f:
		for line in r_f:
			line=line.strip()
			if num_param_per_line == 1:
				return line.split(',')[0]
			else:
				print("Code not written yet for num_param_per_line > 1 !!")
				return None
# ----------------------------------------------------------------------------------------------------------#
def check_inferences(inferDir,requiredTrees,workingDir):

	failedfile = workingDir+'/Status.txt'
	with open(inferDir+'/PassFail_Trees.txt','w') as passFail_f:

		completedPhylogenies=[]
		badPhylogenies=[]
		for tree_dir in requiredTrees:
			completedFile = tree_dir + '/inference_step_complete'
			noDuplFile = tree_dir + '/NO_DUPL'
			noSimulationFile = tree_dir + '/NO_SIM'

			treeNum = re.search("(\d+)$", tree_dir).group(1)
			if os.path.exists(completedFile):
				completedPhylogenies.append(treeNum)
				passFail_f.write("Passed Tree at: %s\n" %tree_dir)
			elif os.path.exists(noDuplFile):
				completedPhylogenies.append(treeNum)
				passFail_f.write("Passed Tree (NoDuple) at: %s\n" %tree_dir)
			elif os.path.exists(noSimulationFile):
				badPhylogenies.append(treeNum)
				passFail_f.write("Bad Tree (no Sim) at: %s\n" % tree_dir)
			else:
				badPhylogenies.append(treeNum)
				passFail_f.write("Failed Tree at: %s\n" % tree_dir)

	passFail_f.close()
	print("check_inferences: completedPhylogenies,badPhylogenies,requiredTrees:")
	print(len(completedPhylogenies))
	print(len(badPhylogenies))
	print(len(requiredTrees))

	#Check if all trees have results:
	#if(len(completedPhylogenies) + len(badPhylogenies) == len(requiredTrees)):
	return completedPhylogenies
	#else:
		# Error accured and 1 or more trees didn't finish the run
	#	return 0

# ----------------------------------------------------------------------------------------------------------#
def summarize_models(working_tree_dirs,modelsSummaryFile):

	with open(modelsSummaryFile,'w') as model_sum_f:
		model_sum_f.write('### Models summary ###\n')
		model_sum_f.write('# Best model per tree:\n')
		models = dict()
		for tree_dir in working_tree_dirs:
			tree_num = re.search('ChromEvol_Tree_(\d+)$', tree_dir).group(1)
			treeBestModel = read_val_from_file(tree_dir + '/ChosenModel_LowAIC.txt', 1)
			model_sum_f.write('tree %s: %s\n' %(tree_num,treeBestModel))
			if treeBestModel in models.keys():
				models[treeBestModel]+=1
			else:
				models[treeBestModel]=1
		model_sum_f.write('\n# Models frequencies:\n')
		for model in models.keys():
			model_sum_f.write('%s: %d\n' %(model,models[model]))


		if 'CONST_RATE_NO_DUPL' not in models.keys():
			models['CONST_RATE_NO_DUPL'] = 0
		if 'LINEAR_RATE_NO_DUPL' not in models.keys():
			models['LINEAR_RATE_NO_DUPL'] = 0
		noDuplModels = models['CONST_RATE_NO_DUPL'] + models['LINEAR_RATE_NO_DUPL']
		if noDuplModels/len(working_tree_dirs) > 0.5 :
			model_sum_f.write('\n# More than 50% NO_DUPL models')
			model_sum_f.close()
			return 0
		else:
			model_sum_f.close()
			return 1
			
# ----------------------------------------------------------------------------------------------------------#
def computeReliabilityProfile(INF, SIM, type):
	# sub compute reliability profile
	# compares simulated ploidy against ploidy inferred from simulations
	REL = {"TP":0, "TN":0, "FP":0, "FN":0}
	# PP threshold
	if type == 1:
		for key in INF.keys():
			if (INF[key] == SIM[key] and INF[key] == 1):
				REL["TP"] +=1
			elif (INF[key] == SIM[key] and INF[key] == 0):
				REL["TN"] +=1
			elif (INF[key] != SIM[key] and INF[key] == 1):
				REL["FP"] +=1
			elif (INF[key] != SIM[key] and INF[key] == 0):
				REL["FN"] +=1
	elif type == 0:
		for key in INF.keys():
			if (INF[key] == SIM[key] and INF[key] == 0):
				REL["TP"] +=1
			elif (INF[key] == SIM[key] and INF[key] == 1):
				REL["TN"] +=1
			elif (INF[key] != SIM[key] and INF[key] == 0):
				REL["FP"] +=1
			elif (INF[key] != SIM[key] and INF[key] == 1):
				REL["FN"] +=1

	return REL
# ----------------------------------------------------------------------------------------------------------#
def parseExpectationsFile(file,thres,ploidyCallType):
	# sub parse expectations file
	# returns expected number of events from root to leaf.
	# Can either work on an inferred expectations file or simulation events file

	# LEAF    GAIN    LOSS    DUPLICATION     DEMI-DUPLICATION	(BASE_NUMBER)
	# Onychostoma_gerlachi    0       1.01103 0       0
	# Tor_macrolepis  0       0.0509317       1.84544 0.000364298
	# Tor_putitora    0       0.0111298       1.99996 6.806e-07

	with open(file,'r') as expFile_f:
		perLeaf_EXPECTATIONS_lines=[]
		valid_lines=0
		for line in expFile_f:
			if "LEAF	GAIN	LOSS	DUPLICATION	DEMI-DUPLICATION" in line:
				valid_lines = 1
				continue
			if valid_lines == 1:
				perLeaf_EXPECTATIONS_lines.append(line.strip())
	expFile_f.close()

	EXP = dict()  # store expectations indexed by species name
	requiredEvents = ploidyCallType.split(' ')
	basenumEvents = None
	for line in perLeaf_EXPECTATIONS_lines:
		split_line = line.split('\t')
		if split_line[3]: duplEvents=float(split_line[3])
		if split_line[4]: demiduplEvents=float(split_line[4])
		if len(split_line) > 5:
			basenumEvents=float(split_line[5])
		# calculate number of events
		events=0
		for callType in requiredEvents:
			if duplEvents and callType == 'DUPL':
				events += duplEvents
			if demiduplEvents and callType == 'DEMI':
				events += demiduplEvents
			if basenumEvents and callType == 'BASE':
				events += basenumEvents
		if events >= thres:
			EXP[split_line[0]] = 1
		else:
			EXP[split_line[0]] = 0

	return EXP
# ----------------------------------------------------------------------------------------------------------#
# sub test threshold
# given a threshold, the function returns a hash with the total
# numbers of TP,TN,FP,FN in all simulations
def test_threshold_for_power(working_tree_dirs, ploidyCallType, thres, type):

	testResults = {"TP":0, "TN":0, "FP":0, "FN":0}
	for tree_dir in working_tree_dirs:
		completedFile = tree_dir + '/inference_step_complete'
		if not os.path.exists(completedFile):
			no_inf_f = open(tree_dir + '/NO_inference_step_complete','w')
			no_inf_f.close()
			continue
		simInferDir = tree_dir + '/simulation/sim_infer/'
		model_dir = [name for name in os.listdir(simInferDir) if os.path.isdir(simInferDir + name)]
		INF_FILE = simInferDir + model_dir[0] + '/expectations.txt'
		# find sim events file
		simulationDir = tree_dir + '/simulation/'
		sim_dirs = [name for name in os.listdir(simulationDir) if os.path.isdir(simulationDir + name)]
		for sim_dir in sim_dirs:
			if 'dir_sim' in sim_dir:
				SIM_FILE = tree_dir + '/simulation/' + sim_dir +'/0/simEvents.txt'
		if not os.path.exists(INF_FILE):
			no_inf_f = open(tree_dir + '/NO_expectations', 'w')
			no_inf_f.write('Missing File: %s' %INF_FILE)
			no_inf_f.close()
			continue
		INF_EXP = parseExpectationsFile(INF_FILE, thres, ploidyCallType)
		SIM_EXP = parseExpectationsFile(SIM_FILE, 1, ploidyCallType)
		result = computeReliabilityProfile(INF_EXP, SIM_EXP, type)
		### get reliability profile
		testResults["TP"] += result["TP"]
		testResults["TN"] += result["TN"]
		testResults["FP"] += result["FP"]
		testResults["FN"] += result["FN"]

	return testResults


# ----------------------------------------------------------------------------------------------------------#
def conf_interval(sortedArr, interval):
	# sub confidence interval
	# recieves an array (reference) and an interval (between 0 and 1) and returns the corresponding percentile
	n = len(sortedArr)
	index = round(n * interval)
	if index > n - 1:
		index -= 1
	bound = sortedArr[index]
	return bound
# ----------------------------------------------------------------------------------------------------------#
def largest_value_hash_median(hash):

	big=0
	best=[]
	for key in hash.keys():
		if hash[key] > big:
			big = hash[key]
			best.append(key)
		elif hash[key] == big:
			best.append(key)
	best.sort()
	median = conf_interval(best,0.5)

	return median
# ----------------------------------------------------------------------------------------------------------#
def determine_PP_threshold_for_power(working_tree_dirs,outFile,ploidyCallType):

	try:
		with open(outFile,'w') as thresfh:
			thresMCC = dict()
			for thresP in np.arange(0.0, 1.01, 0.01):
				thresResultsRef = test_threshold_for_power(working_tree_dirs, ploidyCallType, thresP, 1)
				thresTP=thresResultsRef["TP"]
				thresTN=thresResultsRef["TN"]
				thresFP=thresResultsRef["FP"]
				thresFN=thresResultsRef["FN"]
				mcc = MCC(thresTP, thresTN, thresFP, thresFN)
				thresMCC[thresP] = mcc
				thresfh.write("%f\t%f\n" %(thresP,mcc))
			bestThresP = largest_value_hash_median(thresMCC)
			thresfh.write("#Best threshold: %s\n" %str(bestThresP))
			thresfh.close()
			return bestThresP

	except IOError:
		print("Can't open output thresholds file %s\n" % outFile)
		return

# ----------------------------------------------------------------------------------------------------------#
def determine_DP_threshold_for_power(working_tree_dirs,outFile,ploidyCallType,thresPP):

	try:
		with open(outFile,'w') as thresfh:
			thresMCC = dict()
			for thresD in np.arange(0.0, thresPP, 0.01):
				thresResultsRef = test_threshold_for_power(working_tree_dirs, ploidyCallType, thresD, 1)
				thresTP=thresResultsRef["TP"]
				thresTN=thresResultsRef["TN"]
				thresFP=thresResultsRef["FP"]
				thresFN=thresResultsRef["FN"]
				mcc = MCC(thresTP, thresTN, thresFP, thresFN)
				thresMCC[thresD] = mcc
				thresfh.write("%f\t%f\n" %(thresD,mcc))
			bestthresD = largest_value_hash_median(thresMCC)
			thresfh.write("#Best threshold: %s\n" %str(bestthresD))
			thresfh.close()
			return bestthresD

	except IOError:
		print("Can't open output thresholds file %s\n" % outFile)
		return

# ----------------------------------------------------------------------------------------------------------#
def determine_thresholds_for_power(working_tree_dirs,outFilePP, outFileDP, ploidyCallType):

	## determine PP threshold
	thresPP = determine_PP_threshold_for_power(working_tree_dirs,outFilePP, ploidyCallType)
	## determine DP threshold
	thresDP = determine_DP_threshold_for_power(working_tree_dirs,outFileDP, ploidyCallType, thresPP)

	return thresPP,thresDP

# ----------------------------------------------------------------------------------------------------------#
def parseExpectationsFile_2thres(file,thresPP,thresDP,ploidyCallType):

	# sub parse expectations file with two thresholds
	# similar to parse expectations file, but considers a different threshold for DP and PP
	EXP=dict()

	with open(file,'r') as expFile_f:
		perLeaf_EXPECTATIONS_lines=[]
		valid_lines=0
		for line in expFile_f:
			if "LEAF	GAIN	LOSS	DUPLICATION	DEMI-DUPLICATION" in line:
				valid_lines = 1
				continue
			if valid_lines == 1:
				perLeaf_EXPECTATIONS_lines.append(line.strip())
	expFile_f.close()

	EXP = dict()  # store expectations indexed by species name
	requiredEvents = ploidyCallType.split(' ')
	basenumEvents = None
	for line in perLeaf_EXPECTATIONS_lines:
		split_line = line.split('\t')
		if split_line[3]: duplEvents=float(split_line[3])
		if split_line[4]: demiduplEvents=float(split_line[4])
		if len(split_line)>5:
			basenumEvents=float(split_line[5])
		# calculate number of events
		events=0
		for callType in requiredEvents:
			if duplEvents and callType == 'DUPL':
				events += duplEvents
			if demiduplEvents and callType == 'DEMI':
				events += demiduplEvents
			if basenumEvents and callType == 'BASE':
				events += basenumEvents
		if events >= thresPP:
			EXP[split_line[0]] = 1
		elif events > thresDP and events < thresPP:
			EXP[split_line[0]] = 'NA'
		elif events <= thresDP:
			EXP[split_line[0]] = 0

	return EXP

# ----------------------------------------------------------------------------------------------------------#
def compute_reliability_from_real_for_power(working_tree_dirs, thresPP, thresDP, outFile, ploidyCallType):

	# Run on infers, parse expectations file and add to main hash
	taxaHash=dict()
	for tree_dir in working_tree_dirs:
		resSum = tree_dir + '/result_sum.csv'
		bestModel = read_val_from_file(tree_dir + '/ChosenModel_LowAIC.txt', 1)
		expFile = tree_dir + '/' + bestModel + '/expectations.txt'
		expHash = parseExpectationsFile_2thres(expFile,thresPP,thresDP,ploidyCallType)
		for key in expHash.keys():
			expArr=[]
			if key in taxaHash.keys():
				expArr = taxaHash[key]
				expArr.append(expHash[key])
			else:
				expArr.append(expHash[key])
			taxaHash[key] = expArr

	# Check reliability for each taxon and print to file
	with open(outFile,'w') as OUT:
		taxa = list(taxaHash.keys())
		for taxon in taxa:
			mostFreqKey, majority =majority_in_arr(taxaHash[taxon])
			OUT.write("%s\t%s\t%s\n" %(taxon,mostFreqKey,majority))  # if reliable, print ploidy
	OUT.close()
	return 'DONE'
# ----------------------------------------------------------------------------------------------------------#
def parse_reliability_file(relFile):

	# sub parse reliability file
	print("parse_reliability_file")
	resHash = dict()
	with open(relFile,'r') as relLines_f:
		print(relFile)
		for line in relLines_f:
			line_str = line.strip()
			splitLine = line_str.split('\t')
			info = [splitLine[1],splitLine[2]]
			resHash[splitLine[0]] = info

	relLines_f.close()
	return resHash
# ----------------------------------------------------------------------------------------------------------#
def computeFinalReliabilityProfile(INF,SIM):

	# sub compute final reliability profile
	# compares simulated ploidy against ploidy inferred from simulations
	REL=dict()
	for key in INF.keys():
		if INF[key] == SIM[key]:
			if key in REL.keys():
				REL[key]+=1
			else:
				REL[key]=1
		else:
			if key not in REL.keys():
				REL[key]=0

	return REL
# ----------------------------------------------------------------------------------------------------------#
def compute_reliability_from_sim_for_power(working_tree_dirs,thresDP, thresPP, ploidyCallType,outFile,inferExpRef):

	TOTAL_REL=dict()
	total=0
	for tree_dir in working_tree_dirs:
		# check if can work with this phylogeny
		completedFile = tree_dir + '/inference_step_complete'
		if not os.path.exists(completedFile):
			no_inf_f = open(tree_dir + '/NO_inference_step_complete','w')
			no_inf_f.close()
			continue

		simulationDir = tree_dir + '/simulation/'
		sim_dirs = [name for name in os.listdir(simulationDir) if os.path.isdir(simulationDir + name)]
		for sim_dir in sim_dirs:
			if 'dir_sim' in sim_dir:
				SIM_FILE = tree_dir + '/simulation/' + sim_dir + '/0/simEvents.txt'
		# find sim infer expectations
		simInferDir = tree_dir + '/simulation/sim_infer/'
		model_dir = [name for name in os.listdir(simInferDir) if os.path.isdir(simInferDir + name)]
		INF_FILE = simInferDir + model_dir[0] + '/expectations.txt'
		total+=1

		INF_EXP = parseExpectationsFile_2thres(INF_FILE, thresDP, thresPP, ploidyCallType)
		SIM_EXP = parseExpectationsFile(SIM_FILE, 1, ploidyCallType)
		REL_PROFILE = computeFinalReliabilityProfile(INF_EXP, SIM_EXP)
		for key in REL_PROFILE.keys():
			if key in TOTAL_REL.keys():
				TOTAL_REL[key]+=REL_PROFILE[key]
			else:
				TOTAL_REL[key] = REL_PROFILE[key]

		inferExp = inferExpRef
		with open(outFile,'w') as out_f:
			for key in inferExp.keys():
				relPercent = TOTAL_REL[key]/total
				species=key
				relValues = inferExp[key]
				ploidy = relValues[0]
				out_f.write('%s\t%s\t%s\n' %(str(species),str(ploidy),str(relPercent)))
		out_f.close()

	return 'DONE'
# ----------------------------------------------------------------------------------------------------------#
def parse_counts_file(countsFile):

	counts_dict_all=dict()
	spc_in_couns_list=[]
	flag_name = 0
	with open(countsFile,'r') as counts_f:
		for line in counts_f:
			if not line.strip():  continue   #To skip empty lines
			if '>' in line and flag_name == 0:
				flag_name = 1
				spc_name = line.strip().replace('>','')
				spc_in_couns_list.append(spc_name)
			else:
				if 'x' in line.strip():
					counts_dict_all[spc_name]=line.strip()
				else:
					counts_dict_all[spc_name]=int(line.strip())
				flag_name=0
	return counts_dict_all
# ----------------------------------------------------------------------------------------------------------#
def summarize_reliability_for_power(relSim,relInf,countsFile,outFile):

	# parse files
	relSimHash = parse_reliability_file(relSim)
	relInfHash = parse_reliability_file(relInf)

	# compare and print to output
	with open(outFile,'w') as OUT:
		#OUT.write('"Taxon", "Chromosome count", "Ploidy inference", "Phylogeny robustness score", '
		#'"Simulation reliability score","Combined score"\n')
		OUT.write('Taxon,Chromosome count,Ploidy inference,Phylogeny robustness score,'
		'Simulation reliability score,Combined score\n')

		countsHash = parse_counts_file(countsFile)
		taxa =  list(relSimHash.keys())
		for taxon in taxa:
			count = countsHash[taxon]
			infer = relInfHash[taxon][0]
			phyloRobust = relInfHash[taxon][1]
			simReliability = relSimHash[taxon][1]
			combinedScore = (float(phyloRobust)/2 + float(simReliability)/2)
			OUT.write('"%s","%s","%s","%s","%s","%s"\n' % (taxon,count,infer,phyloRobust,simReliability,combinedScore))

	return
# ----------------------------------------------------------------------------------------------------------#
def getIndexByModel(data, model):
	counter = 0;
	for d in data: 
		if d['name'] == model:
			return counter
		counter += 1
	return -1
# ----------------------------------------------------------------------------------------------------------#
def fail(error_msg, type, error_file_path):
    with open(error_file_path, 'w') as error_f:
        error_f.write(error_msg + '\n')
    raise Exception(error_msg, type)