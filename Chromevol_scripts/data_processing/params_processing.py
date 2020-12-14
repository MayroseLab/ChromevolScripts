from MA_defs import *
from data_processing import best_model


if sanity==1:
	d["mainType"] = "Optimize_Model"
	d["_dataFile"] = data_file
	d["_maxChrNum"] = -10
	d["_minChrNum"] = -1
else:
	d["mainType"] = "mainSimulate"
	d["_freqFile"] = freq_file
	d["_simulationsIter"] = 100
	d["_simulationsJumpsStats"] = "expStats.txt"
	d["_maxChrNumForSimulations"] = max(counts) * 10
	d["_simulationsTreeLength"] = 4


d["_outDir"] = outDir
d["_treeFile"] = tree_file
d["_logFile"] = "log.txt"
d["_logValue"] = 6
d["_maxChrNum"] = -10
d["_minChrNum"] = -1
d["_maxOptimizationIterations"] = 5
d["_epsilonLLimprovement"] = 0.01
d["_optimizeIterNum"] = "0,1,3"
d["_optimizePointsNum"] = "5,2,1"
d["_branchMul"] = 1


def initialize_defaults():
	'''
		initialize parameters default values
	:return: parameters dictionary with fixed parameters
	'''
	d = {}

	d["_maxChrNum"] = -10
	d["_minChrNum"] = -1
	d["_branchMul"] = 999
	d["_simulationsNum"] = 1000
	d["_logFile"] = "log.txt"
	d["_logValue"] = 6
	d["_maxOptimizationIterations"] = 5
	d["_epsilonLLimprovement"] = 0.01
	d["_optimizeIterNum"] = "0,1,3"
	d["_optimizePointsNum"] = "5,2,1"
	d["_simulationsIter"] = 100
	d["_simulationsTreeLength"] = 4

	return d

def create_user_param_dict(filename):
	d = {}
	with open(filename, "r") as params_file:
		for line in params_file:
			line = line.strip()
			name = re.search("(.*)\s(.*)",line).group(1)
			val = re.search("(.*)\s(.*)", line).group(2)
			d[name] = val
	return d

def create_params_dict(outDir, dataFile, treeFile, params_from_user):

	d = initialize_defaults()

	d["_mainType"] = "mainSimulate"
	d["_outDir"] = outDir
	if os.path.isfile(dataFile):
		d["_dataFile"] = dataFile
	d["treeFile"] = treeFile

	if os.path.isfile(params_from_user):
		tmp_d = create_user_param_dict(params_from_user)
	else:
		#########################>>>>>>>>>>>>>>>>>>>> need model for path names, but user don't necessarily specify model >>>>>>>>>>>>>>>######################
		tmp_d = best_model.get_params(main_res_dir + model + CE_res_filename,main_res_dir + model + root_freq_filename,max(counts) * 10) # parse from existing CE results file
		bntpv_vec = create_bntpv(main_res_dir + model_name + expectation_file, main_res_dir + model_name + mlAncTree,d["_baseNumber"])
		d["_baseTransitionProbs"]  = bntpv_vec
		d["_maxChrNumForSimulations"] = max(counts) * 10
	# d = {**d, **tmp_d}
	return d

