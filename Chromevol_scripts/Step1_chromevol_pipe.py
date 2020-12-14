from Bio import Phylo


#_treesFile /groups/itay_mayrose/michaldrori/MSA_JmodelTree/output_MD/Leptinella/Leptinella_Chromevol_prune/parsemb_trees.tre
#_dataFile /groups/itay_mayrose/michaldrori/MSA_JmodelTree/output_MD/Leptinella/Leptinella_Chromevol_prune/Leptinella.counts
#_outDir /groups/itay_mayrose/michaldrori/MSA_JmodelTree/output_MD/Leptinella/Leptinella_Chromevol_prune/chromevol_out
#_name Leptinella
#_paramTemplates /groups/itay_mayrose/michaldrori/Pipeline_Files/power_PARAM_templates/
#_chromevolExe /groups/itay_mayrose/michaldrori/scripts/chromEvol_v2.exe
#_cpusNum 1
#_runModels CONST

#----------------------------------------------------------------------------------------------------------#
def print_to_log(lineData, Log_f):

    with open (Log_f,'w+') as f_log:
        f_log.write(lineData)
    f_log.close()
    return
#----------------------------------------------------------------------------------------------------------#
def read_pip_control(working_dir,pip_file):

    param_dict = {'_topologyNum': 100, '_simNum':100, '_baseNum':0, '_ploidyCallType':'DUPL DEMI BASE',
                           '_cpusNum':0, '_ploidyModelsOnly':1, '_runModels':'ALL', '_excludeModels':0,
                           '_plot':0, '_fixRoot':0, '_name':"chromEvol"}

    with open(pip_file,'r') as pip_f:
        for line in pip_f:
            if '#_dataFile' in line:
                countsFile = line.strip().split(' ')[1]
                param_dict['_dataFile']=countsFile
            if '#_treesFile' in line:
                inTreesFile = line.strip().split(' ')[1]
                param_dict['_treesFile'] = inTreesFile
            if '#_topologyNum' in line:
                infTreeNum = line.strip().split(' ')[1]
                param_dict['_topologyNum'] = infTreeNum
            if '#_simNum' in line:
                simNum = line.strip().split(' ')[1]
                param_dict['_simNum'] = simNum
            if '#_paramTemplates' in line:
                paramTemplateDir = line.strip().split(' ')[1]
                param_dict['_paramTemplates'] = paramTemplateDir
            if '#_baseNum' in line:
                baseNum = line.strip().split(' ')[1]
                param_dict['_baseNum'] = baseNum
            if '#_outDir' in line:
                workdir = line.strip().split(' ')[1]
                param_dict['_outDir'] = workdir
                logFile = workdir + "/log.txt"
            if '#_chromevolExe' in line:
                chromEvolExe = line.strip().split(' ')[1]
                param_dict['_chromevolExe'] = chromEvolExe
            if '#_ploidyCallType' in line:
                ploidyCallType = line.strip().split(' ')[1]
                param_dict['_ploidyCallType'] = ploidyCallType
            if '#_cpusNum' in line:
                cpusNum = line.strip().split(' ')[1]
                param_dict['_cpusNum'] = cpusNum
            if '#_ploidyModelsOnly' in line:
                duplModelsOnly = line.strip().split(' ')[1]
                param_dict['_ploidyModelsOnly'] = duplModelsOnly
            if '#_logFile' in line:
                logFile = line.strip().split(' ')[1]
                param_dict['_logFile'] = logFile
            if '#_runModels' in line:
                runModels = line.strip().split(' ')[1]
                param_dict['_runModels'] = runModels
            if '#_excludeModels' in line:
                excludeModels = line.strip().split(' ')[1]
                param_dict['_excludeModels'] = excludeModels
            if '#_plot' in line:
                plot = line.strip().split(' ')[1]
                param_dict['_plot'] = plot
            if '#_fixRoot' in line:
                fixRoot = line.strip().split(' ')[1]
                param_dict['_fixRoot'] = fixRoot
            if '#_name' in line:
                analysisName = line.strip().split(' ')[1]
                param_dict['_name'] = analysisName

    for param in param_dict.keys():
        line_to_print = ('%s:%s' %(param,param_dict[param]))
        print_to_log(line_to_print,working_dir+'/Log.txt')

    return
#----------------------------------------------------------------------------------------------------------#
def edit_counts_file(consenTree):#, countsFile, outCounts):

    tree = Phylo.read(consenTree, "newick")
    for leaf in tree.get_terminals():
        print(leaf.name)

    return


consenTree ='/groups/itay_mayrose/michaldrori/MD_ChromEvol/SERVER_chromEvol/TreeFile'
edit_counts_file(consenTree)
