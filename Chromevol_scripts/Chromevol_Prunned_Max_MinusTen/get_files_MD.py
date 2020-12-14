"""
	Collects required files for chromEvol run
"""

import os
import sys
import shutil
#from ete3 import Tree
import subprocess

#COUNTS_DIR = '/groups/itay_mayrose/share/ploidb/ccdb/pergenus-6SEP15/'
COUNTS_DIR = '/groups/itay_mayrose/michaldrori/Pipeline_Files/Counts_perGenus/'

#POWER_PARAM_DIR = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/power_PARAM_templates/'
POWER_PARAM_DIR = '/groups/itay_mayrose/michaldrori/Pipeline_Files/power_PARAM_templates/'

#CHROMEVOL_EXE = '/groups/itay_mayrose/itaymay/pupkoSVN/trunk/programs/chromEvol/chromEvol_source-2.0_jack/chromEvol.exe'
CHROMEVOL_EXE = '/groups/itay_mayrose/michaldrori/scripts/chromEvol.exe'


MSA_OUTPUT_PATH = "/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output_MD/"
STEP1_PATH = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/step1_preparation.pl'
STEP1_LILACH_PATH = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/step1_preparation_queue_lilach.pl'


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Run external process
#
def exec_external_command_redirect_output(command_to_exec, log_file, outfile=None, errfile=None, cwd=None):
    out_handle = None
    err_handle = None
    if outfile is None:
        command_out = subprocess.PIPE
    else:
        out_handle = open(outfile, "w")
        command_out = out_handle

    if errfile is None:
        command_err = subprocess.PIPE
    else:
        err_handle = open(errfile, "w")
        command_err = err_handle

    log_file.write("Executing the following command %s - out written to %s error written to %s" % (command_to_exec, outfile, errfile))
    log_file.write("cwd=%s" % cwd)
    p = subprocess.Popen(command_to_exec, shell=True, cwd=cwd, stdout=command_out, stderr=command_err)
    stdout, stderr = p.communicate()
    print("After subprocess, Before wait")
    retval = p.wait()
    print("After wait %s" % retval)

    if retval == 0:
        log_file.write("Execution return code was %i for command %s" % (retval, command_to_exec))
    else:
        log_file.write("Execution return code was %i for command %s" % (retval, command_to_exec))

    if out_handle is not None:
        out_handle.close()
    if err_handle is not None:
        err_handle.close()

    return retval
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def prepare_sh_file(genus,chromEvol_dir,queue):
	genus_dir = chromEvol_dir
	sh_file = genus_dir + 'send_to_queue.sh'
	print(sh_file)
	job_name = genus + '_chromEvol'
	PIP_control = genus_dir + 'PIP_control'
	out_file = genus_dir + 'queue.OU'
	err_file = genus_dir + 'queue.ER'
	if queue == 'itaym':
		script = STEP1_PATH
	elif queue == 'lilach':
		script = STEP1_LILACH_PATH
	text = ('#!/bin/tcsh\n'
			'#$ -N ' + job_name + '\n'
			'#$ -S /bin/tcsh\n'
			'#$ -cwd\n'
			'#$ -e ' + err_file + '\n'
			'#$ -o ' + out_file + '\n'
			#'\n' + TreeMap_command + '\n'
			'module load  perl/perl-5.16.3\n'
			'perl ' + script + ' ' + PIP_control +
			'\n'
	)
	with open(sh_file, 'w') as f:
		f.write(text)

def send_to_queue(genus,chromEvol_dir,queue):
	genus_dir = chromEvol_dir
	sh_file = genus_dir + 'send_to_queue.sh'
	if queue == 'itaym':
		os.system('qsub -l itaym -p -1 ' + sh_file)
	elif queue == 'lilach':
		os.system('qsub -l lilach -p -3 ' + sh_file)

def clean(genus,chromEvol_dir):
	genus_dir = chromEvol_dir
	chromEvol_out = chromEvol_dir + 'chromevol_out'
	#print("XXXXXX clean " + chromEvol_out)
	if os.path.isdir(chromEvol_out):
		shutil.rmtree(chromEvol_out)
		os.remove(genus_dir + 'queue.ER')
		os.remove(genus_dir + 'queue.OU')
		os.remove(genus_dir + 'send_to_queue.sh')


def create_genus_dir(genus,output_dir,prune_flag):
	"""
	Creates directory for genus. Returns False if directory was not created.
	"""
	if prune_flag==1:
		dir_path = output_dir  + genus + "/" + genus + "_Chromevol_prune/"
	elif prune_flag == 0:
		dir_path = output_dir  + genus + "/" + genus + "_Chromevol/"
	if os.path.exists(dir_path):
		#print('Directory for genus',genus,'already exists. Skipping...')
		#return False
		shutil.rmtree(dir_path)
	try:
		os.makedirs(dir_path)
		return True
	except:
		print("Can't create directory for genus",genus,". Skipping...")
		return False

def get_trees(genus,output_dir):
	"""
	Copies trees to output location. Returns False if copying failed.
	"""
	#files_dir = TREES_DIR + genus + '/rnr/' # this is if we're taking trees AFTER rnr
	files_dir = MSA_OUTPUT_PATH + genus + '/concat/' # this is if we're taking trees BEFORE rnr
	#MAP_tree_file = files_dir + 'parsemb_map_tree.tre'
	trees_file = files_dir + 'parsemb_trees.tre'
	out_dir = output_dir  + genus + "/" + genus + "_Chromevol/"
	try:
		shutil.copy(trees_file,out_dir)
	except:
		print("Can't copy trees for genus",genus)
		return False
	return True


def get_trees_prune(genus,output_dir):
	"""
	Copies trees to output location. Returns False if copying failed.
	"""
	#files_dir = TREES_DIR + genus + '/rnr/' # this is if we're taking trees AFTER rnr
	files_dir = MSA_OUTPUT_PATH + genus + '/concat/' # this is if we're taking trees BEFORE rnr
	#MAP_tree_file = files_dir + 'parsemb_map_tree.tre'
	trees_file = files_dir + 'parsemb_trees.tre'
	out_dir = output_dir  + genus + "/" + genus + "_Chromevol_prune/"
	try:
		shutil.copy(trees_file,out_dir)
	except:
		print("Can't copy trees for genus",genus)
		return False
	return True


def get_counts(genus,output_dir,prune_flag):
	#first_letter = genus[0]
	#counts_file = COUNTS_DIR + first_letter + '/' + genus + '.fasta'
	counts_file = COUNTS_DIR + genus + '.counts'
	if prune_flag == 1:
		out_dir = output_dir  + genus + "/" + genus + "_Chromevol_prune/"
	else:
		out_dir = output_dir  + genus + "/" + genus + "_Chromevol/"
	try:
		shutil.copy(counts_file,out_dir)
	except:
		print("Can't copy counts for genus",genus)
		return False
	return True


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a list of species
# Get all the species and their counts as in CCDB
def get_species_ToRemove(genus,output_dir):
	example = 'Campanula_pollinensis'
	Species_withCounts_list=[]; Species_withCounts_count = 0; species_name='Null'
	species_list_OnTree_list=[];species_list_OnTree_count=0
	Species_ToRemove=[];removed_species_count=0
	#Open file to write species names for removal:
	species_remove_list_file = MSA_OUTPUT_PATH + genus + '/' + genus + '_Chromevol_prune/remove_species_list.txt'
	f_species_to_remove = open(MSA_OUTPUT_PATH + genus + '/' + genus + '_Chromevol_prune/remove_species_list.txt','w')
	# Get all species names with counts (according to counts file based on CCDB):
	counts_file = MSA_OUTPUT_PATH + genus + '/' + genus + '_Chromevol_prune/' +genus + '.counts'
	f= open(counts_file,'r')
	for line in f:
		line=line.strip('\n')
		if '>' in line:
			species_name = line.strip('>')
			genus_in_species = species_name.split('_')[0]
			if genus == genus_in_species:
				key_species_name = species_name
				#species_list.append(species_name)
				#Species_withCounts_count+=1
		else:
			count_val = line
			if genus == genus_in_species:
				Species_withCounts_list.append(key_species_name)
				Species_withCounts_count+=1

	# Check which species on Tree are missing chrom count:
	tree_file = MSA_OUTPUT_PATH + genus + '/concat/parsemb_trees.tre'
	f_tree = open (tree_file,'r')
	for line in f_tree:
		split_line = line.replace('(','*')
		split_line = split_line.replace(':','*')
		split_line = split_line.replace(',','*')
		split_line = split_line.split('*')
		print("Line Tree !!!\n")
		break
	# Save species on tree and the total count of species on Tree:
	for item in split_line:
		if genus in item:
			species_list_OnTree_list.append(item)
			species_list_OnTree_count+=1
			if item not in Species_withCounts_list:
				Species_ToRemove.append(item)
				f_species_to_remove.write(item+'\n')
				removed_species_count+=1
	# Check if qualifies to run chromevol: more than 5 species with counts:
	if (species_list_OnTree_count - removed_species_count) < 5:
		f_lessThanFive = open(MSA_OUTPUT_PATH + genus + '/' + genus + '_Chromevol_prune/LessThanFiveWithCounts.txt','w')
		f_lessThanFive.close()
		return
	#print("Remove species from list:")
	#print(Species_ToRemove)
	f_species_to_remove.close()

	if os.stat(species_remove_list_file).st_size == 0:
		get_trees_prune(genus,output_dir)
		return


	concat_dir=MSA_OUTPUT_PATH + genus + '/concat/'
	outfile = MSA_OUTPUT_PATH + genus + '/' + genus + '_Chromevol_prune/parsemb_trees.tre'
	TreePrune_command = ("/share/apps/R301/bin/R CMD BATCH '--args working_dir=" + '"' + concat_dir + '"' + " trees.file=" + '"' + tree_file + '"' + " remove_species_list=" + '"' + species_remove_list_file + '"' + " outfile=" + '"' + outfile + '"' + "'" + " /groups/itay_mayrose/michaldrori/scripts/Rscripts/prune.trees.R " + concat_dir + "prune.trees.Rout " + "r" + genus)
	log_file = MSA_OUTPUT_PATH + genus + '/' + genus + '_Chromevol_prune/pruneTree_log.txt'
	f_log=open(log_file,'w')
	exec_external_command_redirect_output(TreePrune_command, f_log)
	f_log.close()
#	try:
#		os.system(TreePrune_command)
#	except:
#		print("Can't create pruned tree file for %s",genus)
#		return False
	print("After pruned Tree")
	return
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



def create_PIP_control(genus,output_dir,prune_flag):
	if prune_flag == 1:
		genus_dir = output_dir + genus + "/" + genus + "_Chromevol_prune/"
	else:
		genus_dir = output_dir + genus + "/" + genus + "_Chromevol/"
	control_file = genus_dir + 'PIP_control'
	trees_file = genus_dir + 'parsemb_trees.tre'
	counts_file = genus_dir + genus + '.counts'
	out_dir = genus_dir + 'chromevol_out'
	constant_text = ('_paramTemplates ' + POWER_PARAM_DIR + '\n' + '_chromevolExe ' + CHROMEVOL_EXE + '\n' + '_cpusNum 1\n' '_runModels CONST\n')
	#constant_text = ('_paramTemplates /groups/itay_mayrose/michaldrori/MD_ChromEvol/power_PARAM_templates/\n'
	#'_chromevolExe /groups/itay_mayrose/itaymay/pupkoSVN/trunk/programs/chromEvol/chromEvol_source-2.0_jack/chromEvol\n'
	#'_cpusNum 1\n'
	#'_runModels CONST')

	with open(control_file, 'w') as fh:
		fh.write('_treesFile ' + trees_file + '\n')
		fh.write('_dataFile ' + counts_file + '\n')
		fh.write('_outDir ' + out_dir + '\n')
		fh.write('_name ' + genus + '\n')
		fh.write(constant_text + '\n')

def do_it_all(genus,output_dir,prune_flag):
	print(prune_flag)
	create_genus_dir(genus,output_dir,prune_flag)
	get_counts(genus,output_dir,prune_flag)
	create_PIP_control(genus,output_dir,prune_flag)
	if prune_flag == 1:
		get_species_ToRemove(genus,output_dir)
		#get_trees_prune(genus,output_dir,species_to_keep)
	else:
		get_trees(genus,output_dir)


if __name__ == "__main__":
	Prune_flag=1
	genera_list = sys.argv[1]
	queue = sys.argv[2]
	if (len(sys.argv)) == 4:
		Prune_flag = int(sys.argv[3])

	output_dir = MSA_OUTPUT_PATH
	with open(genera_list) as f:
		genera = f.readlines()
	for genus in genera:
		genus = genus.strip()
		working_dir_r = MSA_OUTPUT_PATH + genus + "/concat/"
		print(genus)
		#if os.path.exists(chromEvol_dir + "/chromevol_out/log.txt"):
		#	print("log.txt exists for " + genus + "Check if needs rerun and erase old files")
		#	continue
		#concat_dir = MSA_OUTPUT_PATH + "/" + genus + "/concat/"
		#TreeMap_command = ("/share/apps/R301/bin/R CMD BATCH '--args working_dir=" + '"' + concat_dir + '"' + " g=" + '"' + genus + '"' + "'" + " /groups/itay_mayrose/michaldrori/scripts/Rscripts/run_mb_pipe.R " + concat_dir + "run_mb_pipe.Rout " + "r" + genus)
		#os.system(TreeMap_command)
		print(Prune_flag)
		do_it_all(genus,output_dir,Prune_flag)
		print(genus + ": Chromevol prep Completed successfully, continue to create ploidy.csv stage")
		#Send to queue:
		if Prune_flag == 1:
			chromEvol_dir=MSA_OUTPUT_PATH + "/" + genus + "/" + genus + "_Chromevol_prune/"
			if os.path.exists(MSA_OUTPUT_PATH + genus + '/' + genus + '_Chromevol_prune/LessThanFiveWithCounts.txt'):
				print("Less Than 5 species with counts, skip genus %s\n" % genus)
				continue
		else:
			chromEvol_dir=MSA_OUTPUT_PATH + "/" + genus + "/" + genus + "_Chromevol/"
		if os.path.exists(chromEvol_dir + "/chromevol_out/ploidy.csv"):
			print("ploidy.csv exists for " + genus)
			continue
		else:
			clean(genus,chromEvol_dir)
			prepare_sh_file(genus,chromEvol_dir,queue)
			send_to_queue(genus,chromEvol_dir,queue)


			#concat_dir = MSA_OUTPUT_PATH + "/" + genus + "/concat/"
			#TreeMap_command = ("/share/apps/R301/bin/R CMD BATCH '--args working_dir=" + '"' + concat_dir + '"' + " g=" + '"' + genus + '"' + "'" + " /groups/itay_mayrose/michaldrori/scripts/Rscripts/run_mb_pipe.R " + concat_dir + "run_mb_pipe.Rout " + "r" + genus)
