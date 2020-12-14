from MA_defs import *
#from ete3 import Tree
from data_processing import best_model
#import regex as re

def initialize_defaults(output_dir,max_for_sim,nsims,tree_full_path):
    d = {}
    d["_mainType"] = "mainSimulate"
    d["_outDir"] = output_dir
    d["_treeFile"] = tree_full_path
    d["_simulationsJumpsStats"] =  "expStats.txt" ########## try with/out
    if nsims > 0:
        d["_simulationsIter"] = nsims
    else:
	    d["_simulationsIter"] = 100
    d["_maxChrNumForSimulations"] = max_for_sim  ########## need to see if this matters in runtime (100 vs. 10000)
    d["_simulationsTreeLength"] = 4 ##### parsed from res file
    d["_branchMul"] = 1

    return d


def parse_params_from_res_file(d,paramFile,freqFile,expFile, mlAncTree,model_name,orig_counts): # paramFile is where the parameters are found (e.g., CE results file)
    parameters_dictionary = best_model.get_params(paramFile,freqFile)
    d["_freqFile"] = freqFile
    if parameters_dictionary.get("_baseNumber",0)!=0: # a base number model
	    #d["_baseTransitionProbs"] = create_bntpv(expFile, mlAncTree, parameters_dictionary["_baseNumber"])
        d["_maxBaseTransition"] = max(max(orig_counts) - min(orig_counts), parameters_dictionary.get("_baseNumber"))
        #d["_maxBaseTransition"] = max(orig_counts) - min(orig_counts)
    if model_name=="DysDupDem*":
        d["_demiPloidyR"] =  parameters_dictionary.get("_duplConstR",0)
    d.update(parameters_dictionary)

    return d


def create_control_file(filename, working_dir,output_dir,max_for_sim,nsims,model_name,tree_full_path,orig_counts):
    # 1. d = initialize (outDir, main_res_dir + model_name + expectation_file, main_res_dir + model_name + mlAncTree)
    # 2. parse res file OR receive parameters from user (d, model_name)
    d = initialize_defaults(output_dir,max_for_sim,nsims,tree_full_path)
    CE_res_filename, expectation_file, mlAncTree, root_freq_filename, sim_control, statistics_names = fixed_vars()
    d = parse_params_from_res_file(d,working_dir + CE_res_filename,working_dir + root_freq_filename,working_dir + expectation_file, working_dir + mlAncTree,model_name,orig_counts)

    with open(filename, "w+") as control_file:
        for key in d.keys():
            control_file.write(key + " " + str(d[key]))
            control_file.write("\n")


def run_MA(main_res_dir,filename,working_dir,output_dir,orig_counts,model_name,nsims,tree_full_path,CE_path): #run_MA(output_dir + sim_control,main_res_dir + model, output_dir,original_counts)
    real_max = max(orig_counts)
    init_max_for_sim = max(real_max, min(real_max*10,200))
    # get max chr allowed
    with open(main_res_dir + "/chromEvol.res","r") as res_file:
        for line in res_file:
            line = line.strip()
            tmp = re.search("max chromosome allowed: (\d+)", line)
            if tmp:
                max_allowed = int(tmp.group(1))
                break

    for mult in range(10):  # x is the factor by which we increase the previous _maxChrNumForSimulations
        max_for_sim = 100 * mult + init_max_for_sim
        if (max_for_sim < max_allowed):
            max_for_sim = max_allowed

        create_control_file(filename, working_dir, output_dir, max_for_sim,nsims,model_name,tree_full_path,orig_counts)
        os.system('"' + CE_path + '" ' + filename)
        cmd = "grep -R --include='simEvents.txt' 'Total number of transitions to max chromosome: [^0]' " + output_dir
        tmp = os.system(cmd)
        f = open(working_dir + "/increasing_max_chr.txt", "w")
        f.write("Iteration number " + str(mult) + ", Max number is currently " + str(max_for_sim))
        f.close()
        if tmp != 0: # did not hit upper bound
            break # no need to keep increasing the max number

    return(max_for_sim)

def create_bntpv(expFile, mlAncTree, baseNum):

    with open(expFile,'r') as expFile_f:
        perNode_EXPECTATIONS_lines=[]
        look_for_end=0
        for line in expFile_f:
            if "#ALL EVENTS EXPECTATIONS PER NODE" in line and look_for_end == 0:
                look_for_end = 1    #stop when you find the next #+++++++++++++++++++++++++++++ line
            if look_for_end == 1 and "#++++++++++++++++++++++++" in line:
                break
            if look_for_end == 1:
                if "#ALL EVENTS EXPECTATIONS PER NODE" not in line and "BASE-NUMBER" not in line:
                    perNode_EXPECTATIONS_lines.append(line.strip())
    expFile_f.close()

    btnodes=dict()
    for item in perNode_EXPECTATIONS_lines:
        item_split = item.split()
        node = item_split[0]
        weight = float(item_split[5])
        if weight > 0.1:
            if weight > 1: weight = 1
            btnodes[node] = weight
    if not btnodes:
        return "%d=1.00"%baseNum
    # if for some reason there are no base transition events
    #if not btnodes.values():
    #    baseNum=1
    #    return baseNum

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
    for key in list(transitions.keys()):
        if (transitions[key]/totalTrans < 0.05):
            del transitions[key]

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