from collections import OrderedDict



results_dir = '/bioseq/data/results/chromEvol/1557734293/'

select_model_file = 'Chromevol_selected_models.txt'
chrom_res = 'chromEvol.res'
tree_paths = dict()
tree_sel_models = dict()
with open(results_dir+select_model_file,'r') as models_f:
    for line in models_f:
        # Selected model (AIC) for Tree# 1: BASE_NUM_DUPL,
        line_split = line.split(':')
        selected_model = line_split[1].replace(',','').strip()
        tree_num = line_split[0].split(' ')[-1]
        tree_path = results_dir + 'ChromEvol_Tree_' + str(tree_num) + '/'
        tree_paths[tree_num] =  tree_path
        tree_sel_models[tree_num] =  selected_model


All_Data_dict=dict()
for tree_num in tree_paths.keys():
    #/bioseq/data/results/chromEvol/1557734293/ChromEvol_Tree_4/LINEAR_RATE
    chrom_num_dict=dict()
    gain_loss_dict=dict()
    with open(tree_paths[tree_num] + tree_sel_models[tree_num] + '/' + chrom_res,'r') as res_f:
        strt_gain_loss=0
        for line in res_f:
            if 'chromosome' in line: #min chromosome in data: 7
                line = line.replace('#','')
                line_splt = line.split(' ')
                if 'allowed' in line:
                    chrom_num_type = line_splt[0] + ' chromosome allowed'
                else:
                    chrom_num_type = line_splt[0] + ' chromosome in data'
                chrom_num_dict[chrom_num_type] = line_splt[-1].strip()
            if '#Final Model parameters' in line:
                strt_gain_loss = 1
            if 'F[' in line:
                strt_gain_loss = 0
            if strt_gain_loss == 1 and '#Final Model parameters' not in line:
                line_splt = line.split('\t')
                param_name = line_splt[0]
                param_val = line_splt[1].strip()
                #gain_loss_list.append(line.strip())
                gain_loss_dict[param_name] = param_val

    tree_Data_dict=dict()
    for param in gain_loss_dict.keys():
        tree_Data_dict[param] = gain_loss_dict[param]
    b = OrderedDict(sorted(chrom_num_dict.items()))
    for key in  b.keys():
        tree_Data_dict[key] = b[key]
    All_Data_dict[int(tree_num)] = tree_Data_dict

param_list_all=[]
#build combined params list:
for tree_i in All_Data_dict.keys():
    for param_key in All_Data_dict[tree_i].keys():
        if param_key not in param_list_all:
            param_list_all.append(param_key)

sort_data = OrderedDict(sorted(All_Data_dict.items()))
with open('/groups/itay_mayrose/michaldrori/RESULT_PRINTOUT.csv','w') as out_f:
    out_f.write('"Header1",')
    idx = 2
    for tree_i in sort_data.keys():
        header_tree = tree_i
        out_f.write('"Header%d",' % idx)
        idx+=1
    out_f.write('\n')

    out_f.write('"",')
    for tree_i in sort_data.keys():
        header_tree = tree_i
        out_f.write('"Tree%-15s",' % header_tree)
    out_f.write('\n')
    for param_i in param_list_all:
        out_f.write('"%-25s",' %param_i)
        for tree_i in sort_data.keys():
            out_f.write('"%-15s",' %sort_data[tree_i][param_i])
        out_f.write('\n')



###
###
###
####Model: GENERAL_CHR_MODEL
####Branching Model: GRADUAL
####Input data file: /bioseq/data/results/chromEvol/1557734293/countsFile_edit
####Input tree file: /bioseq/data/results/chromEvol/1557734293/ChromEvol_Tree_1/tree_1
####Tree rooted at N1	sons of root are: N2, Diarrhena_americana,
####total tree length = 6
####min chromosome in data: 7
####min chromosome allowed: 6
####max chromosome in data: 42
####max chromosome allowed: 52
####Half_duplication rate is zero
####Root frequencies: ROOT_LL
####For optimization issues the tree branches were multiplied by 9.27887
####To preserve the original time unit the model parameters should be multiplied by the same factor !!!
####Final Model parameters
###LOSS_CONST	1.38669e-10
###GAIN_CONST	1.38669e-10
###BASE_NUMBER_R	1.25871
###BASE_NUMBER	7
###F[6]=1e-10   F[7]=1   F[8]=1e-10   F[9]=1e-10   F[10]=1e-10   F[11]=1e-10   F[12]=1e-10   F[13]=1e-10   F[14]=1e-10   F[15]=1e-10   F[16]=1e-10   F[17]=1e-10   F[18]=1e-10   F[19]=1e-10   F[20]=1e-10   F[21]=1e-10   F[22]=1e-10   F[23]=1e-10   F[24]=1e-10   F[25]=1e-10   F[26]=1e-10   F[27]=1e-10   F[28]=1e-10   F[29]=1e-10   F[30]=1e-10   F[31]=1e-10   F[32]=1e-10   F[33]=1e-10   F[34]=1e-10   F[35]=1e-10   F[36]=1e-10   F[37]=1e-10   F[38]=1e-10   F[39]=1e-10   F[40]=1e-10   F[41]=1e-10   F[42]=1e-10   F[43]=1e-10   F[44]=1e-10   F[45]=1e-10   F[46]=1e-10   F[47]=1e-10   F[48]=1e-10   F[49]=1e-10   F[50]=1e-10   F[51]=1e-10   F[52]=1e-10
###LogLikelihood = -121.038
###Likelihood = 2.71492e-53
###AIC (Akaike information criterion) = 250.076
