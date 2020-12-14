
"""
	Compares ploidy inferences to ploidy level data found in Kew C-value database
	Modified by Michal SEP 2015
"""

__author__ = 'ItayMNB5'

import pickle
from prepare_traits_file import get_ploidy_data
import sys
import os
import csv
import sqlite3
import re
import statistics
from collections import defaultdict


### GLOBALS
COUNTS_PATH = '/groups/itay_mayrose/michaldrori/PROJECTS/PloiDB_pipeline/Counts_dir/'
phy_robust_dict = {};simu_reli_dict={};comb_scor_dict={}
DEBUG_FLAG = 0
Trees_NUM = 10

#/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output/Silene/Silene_Chromevol/chromevol_out
#KEW_DATA_CCDB='/groups/itay_mayrose/michaldrori/Pipeline_Files/Kew_Ploidy_Median_CCDB1p40.csv'	# SELECT
OUTPUT_DIR = '/groups/itay_mayrose/michaldrori/PROJECTS/PloiDB_pipeline/OUTPUT/'

def get_chosen_model(genus,chrom_dir):
	model_count_max=0
	model_name_select='No_Model'
	save_model_lines=0
	model_summary_file = OUTPUT_DIR + genus + '/' + chrom_dir + '/models_summary'
	f=open(model_summary_file, 'r')
	for line in f:
		if '# Models frequencies:' in line:
			save_model_lines=1
		elif save_model_lines==1:
			line=line.split(':')
			model_name = line[0]
			print(model_name)
			model_count = line[1].strip('\n')
			if DEBUG_FLAG==1: print(int(model_count))
			if int(model_count) > model_count_max:
				model_name_select = model_name
	if DEBUG_FLAG==1: print(model_name_select)
	return model_name_select


def calc_min_max_values(species_list,ploidy_dict,chromCount_dict):
	MinPolyVal=99999;MinDipVal=99999
	MaxDipVal=0;MaxPolyVal=0
	for species_name in species_list:
		if ploidy_dict[species_name] != 'NA' and chromCount_dict[species_name] != 'x':
			if int(ploidy_dict[species_name]) == 1:
				if int(chromCount_dict[species_name]) < (MinPolyVal):
					MinPolyVal = int(chromCount_dict[species_name])
				if int(chromCount_dict[species_name]) > (MaxPolyVal):
					MaxPolyVal = int(chromCount_dict[species_name])
			elif int(ploidy_dict[species_name]) == 0:
				if int(chromCount_dict[species_name]) > (MaxDipVal):
					MaxDipVal = int(chromCount_dict[species_name])
				if int(chromCount_dict[species_name]) < (MinDipVal):
					MinDipVal = int(chromCount_dict[species_name])
	return MinPolyVal,MaxPolyVal,MaxDipVal,MinDipVal


def get_ploidy_data_all(genus,ploidy_file,file_handle,chrom_dir):
	Species_with_Ploidy=0;Species_with_Counts=0
	ploidy_statist=[]
	species_list = [];
	ploidy_dict = {}
	chromCount_dict = {}

	na_num= 0;diploid_num= 0;ploidy_num = 0
	# No Duple case - no ploidy file:
	if(ploidy_file == 'NONE'):
		edit_counts_file = OUTPUT_DIR + genus +'/'+ genus + chrom_dir + '/countsFile_edit'
		chromCount_dict=dict()
		ploidy_dict=dict()
		f= open(edit_counts_file,'r')
		for line in f:
			line=line.strip('\n')
			if '>' in line:
				species_name = line.strip('>')
				species_name_genus = species_name.split('_')[0]
				if genus == species_name_genus:
					key_species_name = species_name
			else:
				count_val = line
				if genus == species_name_genus:
					species_list.append(key_species_name)
					chromCount_dict[key_species_name]=count_val
					if chromCount_dict[key_species_name] != 'x':
						Species_with_Counts+=1
					ploidy_dict[key_species_name]=0	# All Diploids
					phy_robust_dict[key_species_name] = 'NA'
					simu_reli_dict[key_species_name] = 'NA'
					comb_scor_dict[key_species_name] ='NA'
	else:
		with open(ploidy_file,'r') as f:
			reader = csv.DictReader(f)
			indx=0
			for taxon_line in reader:
				print(taxon_line)
				species_name=taxon_line['Taxon']
				phy_robust_dict[species_name] = taxon_line['Phylogeny robustness score']
				simu_reli_dict[species_name] = taxon_line['Simulation reliability score']
				comb_scor_dict[species_name] = taxon_line['Combined score']
				if genus not in species_name:
					continue
				if DEBUG_FLAG == 1: print(">>>>>>>>>>>>>>>>>>>>>>>>>>>")
				if DEBUG_FLAG == 1: print(species_name)
				species_list.append(species_name)
				chromCount_dict[species_name] = taxon_line['Chromosome count']
				if chromCount_dict[species_name] != 'x':
					Species_with_Counts+=1
				ploidy_dict[species_name] = taxon_line['Ploidy inference']
				if ploidy_dict[species_name] != 'NA':
					Species_with_Ploidy+=1
				indx+=1
	#MinPolyVal,MaxPolyVal,MaxDipVal,MinDipVal = calc_min_max_values(species_list,ploidy_dict,chromCount_dict)
	return (species_list,ploidy_dict,chromCount_dict,phy_robust_dict,simu_reli_dict,comb_scor_dict)


# Get all the species and their counts as in CCDB
def get_all_counts(genus,counts_file):
	Species_Count=0
	species_list = []
	chromCount_dict=dict()
	ploidy_dict=dict()
	f= open(counts_file,'r')
	for line in f:
		line=line.strip('\n')
		if '>' in line:
			species_name = line.strip('>')
			genus_in_species = species_name.split('_')[0]
			if genus == genus_in_species:
				key_species_name = species_name
				species_list.append(species_name)
				Species_Count+=1
		else:
			count_val = line
			if genus == genus_in_species:
				chromCount_dict[key_species_name]=count_val
	return species_list,Species_Count,chromCount_dict


def get_BaseNum(genus,chrom_dir):
	model = get_chosen_model(genus,chrom_dir)
	tree_num=1; base_num_total=[]
	while tree_num <= Trees_NUM:
		#"/groups/itay_mayrose/michaldrori/PROJECTS/PloiDB_pipeline/OUTPUT/Acmella/Chromevol_dir/ChromEvol_Tree_2/BASE_NUM"
		mlAncestors_file = OUTPUT_DIR + genus + '/' + chrom_dir + '/ChromEvol_Tree_' + str(tree_num) + '/' + model + '/mlAncestors.tree'
		tree_num+=1
		f= open(mlAncestors_file,'r')
		for line in f:
			base_num = line.split('[N1-', 1)[1].split(']')[0]
		base_num_total.append(int(base_num))
	base_num_avg = statistics.median(base_num_total)
	return base_num_avg,model


def calc_min_max_counts(MinPolyVal,MaxPolyVal,MinDipVal,MaxDipVal):
	minCount=0;maxCount=0
	if MinPolyVal == 99999 and MinDipVal != 99999:
		minCount=min(MaxPolyVal,MinDipVal,MaxDipVal)
	elif MinPolyVal != 99999 and MinDipVal == 99999:
		minCount=min(MinPolyVal,MaxPolyVal,MaxDipVal)
	if MaxPolyVal == 0 and MaxDipVal != 0:
		minCount=min(MinPolyVal,MinDipVal,MaxDipVal)
	elif MaxPolyVal != 0 and MaxDipVal == 0:
		minCount=min(MinPolyVal,MaxPolyVal,MinDipVal)
	return minCount,maxCount


# Check if species will be clasified as diploid or polyploid:
# According to new Ploidy-Criteria files:
def check_ploidy(MinPolyVal,MaxPolyVal,MaxDipVal,MinDipVal,Specie_chromCount):
	if Specie_chromCount == 'x' or Specie_chromCount == 'X':
		return 'NA', 0
	else:
		Specie_chromCount = int(Specie_chromCount)
	suspectPoly_flag=0;ploidy_inf='NA'
	if DEBUG_FLAG == 1: print("$$$ %d,%d,%d,%d\n" %(MinPolyVal,MaxPolyVal,MaxDipVal,MinDipVal))
	# All NA -> All NA
	if (MinPolyVal==99999 and MaxPolyVal==0 and MaxDipVal==0 and MinDipVal==99999):
		ploidy_inf = 'NA'
	# All Polyploids - > All NA
	elif MaxDipVal==0 and MinDipVal==99999:
		ploidy_inf = 'NA'
	# Only Diploids -> Specie_chromCount <= MaxDipVal -> Diploid
	elif MinPolyVal==99999 and MaxPolyVal==0:
		if Specie_chromCount <= MaxDipVal:
			ploidy_inf = '0'
		# In case there are no polyploids and the count is > 1.5MaxDip, then it will be set as polyploid:
		else:
			if Specie_chromCount > (1.5*MaxDipVal):
				ploidy_inf = '1'
			else:
				ploidy_inf = 'NA'
		#Suspect Poly Flag - in case the count*1.5 <= MinDip: set as diploid but falg is 1
		if Specie_chromCount*1.5 <= MinDipVal:
			suspectPoly_flag=1
	#Both Diploids and polyploids:
	else:
		if Specie_chromCount >= MinPolyVal and Specie_chromCount > MaxDipVal:
			ploidy_inf = '1'
		elif Specie_chromCount < MinPolyVal and Specie_chromCount <= MaxDipVal:
			ploidy_inf = '0'
		#Suspect Poly Flag - in case the count*1.5 <= MinDip: set as diploid but falg is 1
		if Specie_chromCount*1.5 <= MinDipVal:
			suspectPoly_flag=1
	return ploidy_inf,suspectPoly_flag

# Check number of species that began Alignment and the number that are included in the concat file, also number of species with diff genus name
def returnSpeciesNumberSeqs(genusName):
	#Count the number of species in the input allseq file and compare to the number in the alignment:
	#Input file:
	species_withNoGenusName=[]
	species_names=[]
	organs_allseq_file = OUTPUT_DIR + genusName + '/' + genusName + '-allseq.fasta'
	if os.path.exists(organs_allseq_file):
		with open(organs_allseq_file,'r') as f:
			for line in f:
				if '>gi' in line:
					parsed_line = line.split('|')
					species_names.append(parsed_line[5])
					if not genusName in parsed_line[5]:
						species_withNoGenusName.append(parsed_line[5])
	species_count_init = (len(set(species_names)))
	species_withNoGenusName=set(species_withNoGenusName)

	species_names_final=[]
	align_file = OUTPUT_DIR + genusName + '/concat/' + genusName + '-concat-aligned.fasta'
	if os.path.exists(align_file):
		with open(align_file,'r') as f:
			for line in f:
				if '>' in line:
					parsed_line = line.split('>')
					species_names_final.append(parsed_line[1])
					#print(species_names_final)
	species_count_final = (len(set(species_names_final)))
	return (species_count_init, species_count_final, species_withNoGenusName)

#----------------------------------------------------------------------------------------------------------#
def create_dir_if_not_exists(dir):
	if not os.path.exists(dir):
		os.makedirs(dir)
		#print_to_log("Created Directory: %s" %dir, working_dir + '/Log.txt')

#----------------------------------------------------------------------------------------------------------#

# This script should cover all:
#-------------------------------------------------------------------------------------------
# add species that are not on the tree: if the species count is within the range of either the dp or the pp in the genus (but not both),
# infer its ploidy level accordingly.
# table with: genus, accepted.num.sp, num.sp.on.tree, num.sp.counts, pct.ploidy.inference, pct.pp (num.pp/num.ploidy), pct.pp (adding sp.not.on.tree)
#-------------------------------------------------------------------------------------------
# Use the following criteria:
# ->  it has a count and its phylogeny reliability is >= 0.9
# ->  it has no count and its combined reliability is >= 0.9
# ->  assign NA otherwise

if __name__ == "__main__":

	#Path_for_Results = '/groups/itay_mayrose/michaldrori/Pipeline_Files/'
	#Path_for_Results = '/groups/itay_mayrose/michaldrori/For_Anna_Ploidb_dist/'

	Prune_Flag=0
	genera_list_file = sys.argv[1].strip()
	Path_for_Results = sys.argv[2] #'/groups/itay_mayrose/michaldrori/PROJECTS/PloiDB_pipeline/Summary_files/'
	# Which data to work on:
	if Prune_Flag==1:
		chrom_dir = '_Chromevol_prune/'
		SpeciesGeneraPloidy_path = Path_for_Results + '/ChromevolSummary_Pruned/'
	else:
		#chrom_dir = '_Chromevol/'
		chrom_dir = 'Chromevol_dir/'
		SpeciesGeneraPloidy_path = Path_for_Results + '/ChromevolSummary/'
	create_dir_if_not_exists(SpeciesGeneraPloidy_path)
	InTree_Flag = dict()
	#error_file = SpeciesGeneraPloidy_path+'Genus_Species_table/Species_Table_ErrFile.txt' #sys.argv[2]
	#out_file = SpeciesGeneraPloidy_path+'Genus_Species_table/Species_Table.csv' #sys.argv[2]
	error_file = SpeciesGeneraPloidy_path+'Species_Table_ErrFile.txt' #sys.argv[2]
	out_file = SpeciesGeneraPloidy_path+'Species_Table.csv' #sys.argv[2]
	out_perGenus_file = SpeciesGeneraPloidy_path+'Genus_AddedSpeciesAnalysis.csv' #sys.argv[2]
	NoPhylogen_GeneraList_file = SpeciesGeneraPloidy_path+'NoPhylogen_Genera.txt'
	NEW_PLOIDY_PATH = SpeciesGeneraPloidy_path+'GenusPloidy_WithAddedSpecies/'
	#Locate in correct path under genus:
	CRITERIA_PLOIDY_PATH = SpeciesGeneraPloidy_path+'GenusPloidy_WithAddedSpecies_withCriteria/'
	f_out_species = open(out_file,'w')
	f_out_genus = open(out_perGenus_file, 'w')
	f_NoPhylogen_GeneraList = open(NoPhylogen_GeneraList_file,'w')
	# genus, accepted.num.sp, num.sp.on.tree, num.sp.counts, pct.ploidy.inference, pct.pp (num.pp/num.ploidy), pct.pp (adding sp.not.on.tree)
	# Missing Columns: Prct_Ploidy_Inference,Prct_Poly (numPoly/numPloidy),Prct_Ploidy_NotOnTree,
	#f_out.write('Genus,Species_Name,Chrom_Count, Ploidy, On_Tree_Flag, \n')  # Gain/Loss same branch expectations.txt file in chosen model
	f_out_species.write('Genus,Taxon,Chromosome_count,Ploidy_inference,Ploidy_inference_Original,Phylogeny_robustness_score,Simulation_reliability_score,Combined_score,On_Tree_Flag,SuspectPoly_Flag\n')
	f_out_genus.write('Genus,in_tree_species,out_of_tree_species,total_species_count,in_tree_poly,in_tree_dip,in_tree_undefined,out_tree_poly,out_tree_dip,out_tree_undefined,out_tree_suspectPoly \n')
	#Ger Kew data
	#Kew_Ploidy_dict,Kew_Counts_dict,DataNotValid_Species = get_ploidy_counts_fromKew()


	# Get all data from Chromevol:
	#Check if ploidy data or noDuple:
	with open(genera_list_file,'r') as f:
		genera_list = [g.strip() for g in f.readlines()]
	count=0
	for genus in genera_list:
		#Check if less than 5 species with counts:
		chrom_PreRun = '_Chromevol_PreRun'
		PreRun_Dir = OUTPUT_DIR + genus +'/'+ genus + chrom_PreRun    #/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output_MD/Anisacanthus/Anisacanthus_Chromevol_PreRun
		if (os.path.exists(PreRun_Dir+'/LessThanFiveWithCounts.txt') and Prune_Flag==1):
			continue
		chromCount_dict=dict()
		species_list=[]
		print(genus)
		#Check if chromevol dir exists:
		if (not os.path.exists(OUTPUT_DIR + genus +'/'+ chrom_dir +'/')):
			f_err = open(error_file,'w')
			f_err.write(genus + ' - No chromevol_out Dir -> ' + OUTPUT_DIR + genus +'/'+ chrom_dir)
			continue
		file_new_ploidy = NEW_PLOIDY_PATH + genus + '_NewPloidy.csv'
		file_criteria_ploidy = OUTPUT_DIR + genus +'/'+ chrom_dir + '/' + genus + '_CriteriaPloidy.csv'
		f_genus = open(file_new_ploidy, 'w')
		f_criteria_ploidy = open(file_criteria_ploidy, 'w')
		f_genus.write('Taxon,Chromosome count,Ploidy inference,Phylogeny robustness score,Simulation reliability score,Combined score,On_Tree_Flag\n')  # Gain/Loss same branch expectations.txt file in chosen model
		f_criteria_ploidy.write('Taxon,Chromosome count,Ploidy inference,Ploidy inference Original,Phylogeny robustness score,Simulation reliability score,Combined score,On_Tree_Flag,SuspectPoly Flag\n')  # Gain/Loss same branch expectations.txt file in chosen model)
		SpeciesNum_Init, SpeciesNum_final,foreign_Species = returnSpeciesNumberSeqs(genus)

		#/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output/Amomum/Amomum_Chromevol/chromevol_out
		ploidy_file = OUTPUT_DIR + genus +'/'+ chrom_dir + '/ploidy.csv'
		noDupl_file = OUTPUT_DIR + genus +'/'+ chrom_dir + 'chromevol_out/NoDuple_MoreThan50.txt'
		NoPhylogen = OUTPUT_DIR + genus +'/'+ chrom_dir + '/Less_than_50_phylogenies.txt'
		NoPhylogenNoDuple = OUTPUT_DIR + genus +'/'+ chrom_dir + '/Less_than_50_phylogenies_NoDUPL.txt'
		if os.path.exists(NoPhylogen):
			f_NoPhylogen_GeneraList.write(genus+'\n')
			continue
		if os.path.exists(NoPhylogenNoDuple):
			f_NoPhylogen_GeneraList.write(genus+',NoDuple\n')
			continue
		#Case: Ploidy file exists:
		if os.path.exists(ploidy_file):
			species_list,ploidy_dict,chromCount_dict,phy_robust_dict,simu_reli_dict,comb_scor_dict = get_ploidy_data_all(genus,ploidy_file,f_out_species,chrom_dir)
			if DEBUG_FLAG==1: print(ploidy_dict)
		#Case: for more than 90 trees, a No Duple model was selected:
		elif os.path.exists(noDupl_file):
			species_list,ploidy_dict,chromCount_dict,phy_robust_dict,simu_reli_dict,comb_scor_dict = get_ploidy_data_all(genus,'NONE',f_out_species,chrom_dir)
		else:
			f_genus.write("NA,NA,NA,NA,NA,NA,NA,NA\n")
			f_criteria_ploidy.write("NA,NA,NA,NA,NA,NA,NA,NA,NA\n")
			f_out_species.write("%s,NoData\n" %genus)
			continue



		SpeciesCount=len(species_list)
		# Collect all counts for genus from .counts file:
		counts_file = COUNTS_PATH + genus + '_counts'
		species_list_allCounts,Species_Num,chromCount_dict_all = get_all_counts(genus,counts_file)

		#Prepare Criteria Ploidy inference:
		ploidyCriteria_dict={}
		for species_name in species_list:
			chrom_count = str(chromCount_dict[species_name])
			if chrom_count == 'x':
				if (comb_scor_dict[species_name] is 'NA') or (float(comb_scor_dict[species_name]) >= 0.9):
					ploidyCriteria_dict[species_name]=ploidy_dict[species_name]
				else:
					ploidyCriteria_dict[species_name]='NA'
			else:
				if (comb_scor_dict[species_name] is 'NA') or (float(comb_scor_dict[species_name]) >= 0.9):
					ploidyCriteria_dict[species_name]=ploidy_dict[species_name]
				else:
					ploidyCriteria_dict[species_name]='NA'
		MinPolyVal,MaxPolyVal,MaxDipVal,MinDipVal = calc_min_max_values(species_list,ploidyCriteria_dict,chromCount_dict)
		print("Min Max Values: %d,%d,%d,%d" %(MinPolyVal,MaxPolyVal,MaxDipVal,MinDipVal))

		# Compare species lists: to check if species in chromevol output:
		out_of_tree_species=0;out_tree_poly=0;out_tree_dip=0;out_tree_undefined=0;out_tree_suspectPoly=0
		in_tree_species=0;in_tree_poly=0;in_tree_dip=0;in_tree_undefined=0
		index=0
		while index < Species_Num:
			species_name = species_list_allCounts[index]
			#Check if species is missing from chromevol output:
			if species_name not in species_list:
				if DEBUG_FLAG==1: print(species_name)
				ploidy_inf,suspectPoly_flag = check_ploidy(MinPolyVal,MaxPolyVal,MaxDipVal,MinDipVal,(chromCount_dict_all[species_name]))
				#Taxon,Chromosome count,Ploidy inference,Phylogeny robustness score,Simulation reliability score,Combined score,On_Tree_Flag
				print(species_name)
				print(ploidy_inf)
				f_genus.write("%s,%s,%s,NA,NA,NA,0,%d\n" %(species_name,chromCount_dict_all[species_name],ploidy_inf,suspectPoly_flag))
				f_criteria_ploidy.write("%s,%s,%s,NA,NA,NA,NA,0,%d\n" %(species_name,chromCount_dict_all[species_name],ploidy_inf,suspectPoly_flag))
				f_out_species.write("%s,%s,%s,%s,NA,NA,NA,NA,0,%d\n" %(genus,species_name,chromCount_dict_all[species_name],ploidy_inf,suspectPoly_flag))
				if ploidy_inf == '0': out_tree_dip+=1
				elif ploidy_inf == '1': out_tree_poly+=1
				elif ploidy_inf == 'NA': out_tree_undefined+=1
				if suspectPoly_flag == 1: out_tree_suspectPoly+=1
				InTree_Flag[species_name]=0
				out_of_tree_species+=1
			index+=1


		for species_name in species_list:
			#f_out.write('%s,%s,%s,' % (genus,str(SpeciesNum_Init),str(SpeciesNum_final)))
			chrom_count = str(chromCount_dict[species_name])
			if chrom_count == 'x':
				if (comb_scor_dict[species_name] is 'NA') or (float(comb_scor_dict[species_name]) >= 0.9):
					f_criteria_ploidy.write("%s,NA,%s,%s,%s,%s,%s,1,0\n" %(species_name,ploidy_dict[species_name],ploidy_dict[species_name],phy_robust_dict[species_name],simu_reli_dict[species_name],comb_scor_dict[species_name]))
					f_out_species.write("%s,%s,NA,%s,%s,%s,%s,%s,1,0\n" %(genus,species_name,ploidy_dict[species_name],ploidy_dict[species_name],phy_robust_dict[species_name],simu_reli_dict[species_name],comb_scor_dict[species_name]))
				else:
					f_criteria_ploidy.write("%s,NA,NA,%s,%s,%s,%s,1,0\n" %(species_name,ploidy_dict[species_name],phy_robust_dict[species_name],simu_reli_dict[species_name],comb_scor_dict[species_name]))
					f_out_species.write("%s,%s,NA,NA,%s,%s,%s,%s,1,0\n" %(genus,species_name,ploidy_dict[species_name],phy_robust_dict[species_name],simu_reli_dict[species_name],comb_scor_dict[species_name]))
				f_genus.write("%s,NA,%s,%s,%s,%s,1,0\n" %(species_name,ploidy_dict[species_name],phy_robust_dict[species_name],simu_reli_dict[species_name],comb_scor_dict[species_name]))
				#f_out_species.write("%s,%s,NA,%s,1\n" %(genus,species_name,ploidy_dict[species_name]))
			else:
				if (comb_scor_dict[species_name] is 'NA') or (float(comb_scor_dict[species_name]) >= 0.9):
					f_criteria_ploidy.write("%s,%s,%s,%s,%s,%s,%s,1,0\n" %(species_name,chromCount_dict[species_name],ploidy_dict[species_name],ploidy_dict[species_name],phy_robust_dict[species_name],simu_reli_dict[species_name],comb_scor_dict[species_name]))
					f_out_species.write("%s,%s,%s,%s,%s,%s,%s,%s,1,0\n" %(genus,species_name,chromCount_dict[species_name],ploidy_dict[species_name],ploidy_dict[species_name],phy_robust_dict[species_name],simu_reli_dict[species_name],comb_scor_dict[species_name]))
				else:
					f_criteria_ploidy.write("%s,%s,NA,%s,%s,%s,%s,1,0\n" %(species_name,chromCount_dict[species_name],ploidy_dict[species_name],phy_robust_dict[species_name],simu_reli_dict[species_name],comb_scor_dict[species_name]))
					f_out_species.write("%s,%s,%s,NA,%s,%s,%s,%s,1,0\n" %(genus,species_name,chromCount_dict[species_name],ploidy_dict[species_name],phy_robust_dict[species_name],simu_reli_dict[species_name],comb_scor_dict[species_name]))
				f_genus.write("%s,%s,%s,%s,%s,%s,1,0\n" %(species_name,chromCount_dict[species_name],ploidy_dict[species_name],phy_robust_dict[species_name],simu_reli_dict[species_name],comb_scor_dict[species_name]))
				#f_out_species.write("%s,%s,%s,%s,1\n" %(genus,species_name,chromCount_dict[species_name],ploidy_dict[species_name]))
			ploidy_inf = str(ploidy_dict[species_name])
			print(ploidy_dict[species_name])
			if ploidy_inf == '0': in_tree_dip+=1
			elif ploidy_inf == '1': in_tree_poly+=1
			elif ploidy_inf == 'NA' : in_tree_undefined+=1
			InTree_Flag[species_name]=1
			in_tree_species+=1
		f_genus.close()
		f_criteria_ploidy.close()
		total_species_count = in_tree_species+out_of_tree_species
		f_out_genus.write("%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n" %(genus,in_tree_species,out_of_tree_species,total_species_count,in_tree_poly,in_tree_dip,in_tree_undefined,out_tree_poly,out_tree_dip,out_tree_undefined,out_tree_suspectPoly))
		##print(InTree_Flag)
	f_out_species.close()