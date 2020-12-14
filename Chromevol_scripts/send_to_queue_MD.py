"""
	Prepare sh files and send to queue for new chromEvol run
"""

import os
import sys
import shutil
#import PipeLine_Paths

# GLOBALS
#chromEvol_dir = '/groups/itay_mayrose/share/PipeLine/3_ChromEvol/ErrWhileRun_ChroOut/'

MSA_OUTPUT_PATH = '/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output'
STEP1_PATH = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/step1_preparation.pl'
STEP1_LILACH_PATH = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/step1_preparation_queue_lilach.pl'

def prepare_sh_file(genus,chromEvol_dir,queue):
	genus_dir = chromEvol_dir
	sh_file = genus_dir + 'send_to_queue.sh'
	print(sh_file)
	job_name = 'chromEvol_' + genus
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
		os.system('qsub -l itaym -p -2 ' + sh_file)
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

if __name__ == "__main__":
	Prune_flag=0
	genera_list = sys.argv[1]
	#chromEvol_dir = sys.argv[2]
	queue = sys.argv[2]

	if (len(sys.argv)) == 4:
		Prune_flag = int(sys.argv[3])

	#print(genera_list,chromEvol_dir,queue)
	with open(genera_list) as f:
		genera = f.readlines()
	for genus in genera:
		genus = genus.strip()
		print(genus)
		if Prune_flag == 1:
			chromEvol_dir=MSA_OUTPUT_PATH + "/" + genus + "/" + genus + "_Chromevol_prune/"
		else:
			chromEvol_dir=MSA_OUTPUT_PATH + "/" + genus + "/" + genus + "_Chromevol/"
		if os.path.exists(chromEvol_dir + "/chromevol_out/ploidy.csv"):
			print("ploidy.csv exists for " + genus)
			sys.exit(0)
		else:
			clean(genus,chromEvol_dir)
			prepare_sh_file(genus,chromEvol_dir,queue)
			send_to_queue(genus,chromEvol_dir,queue)

