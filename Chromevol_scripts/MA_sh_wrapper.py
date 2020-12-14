import os
import argparse
import csv
import utils

from data_processing import best_model

### ARGS
parser = argparse.ArgumentParser(description="produce sh files for running model adequacy, either in regular or sanity mode")
parser.add_argument('--nsims', '-n', help='Number of simulations',required=True)
parser.add_argument('--sanity', '-s', help='Regular mode = 0, sanity mode = 1',required=True, default = 0)
parser.add_argument('--genera', '-f', help='Genera file',required=True)
parser.add_argument('--script', '-c', help='Script path',required=True)
parser.add_argument('--models_flag', '-m', help='models flag',required=True, default = 0) # options: ALL, BEST, OTHERS, DEFINED
parser.add_argument('--models_defined', '-md', help='define which models to run',required=False) # if DEFINED give model's name(s)
parser.add_argument('--results_flag', '-r', help='results flag',required=True, default = 0) # use previous simulations results

args = parser.parse_args()
nsims = int(args.nsims)
sanity = int(args.sanity)
filename = args.genera
script = args.script
models_flag = args.models_flag
models_defined = args.models_defined
results_flag = int(args.results_flag)

lang = "python"
all_models = ["BASE_NUM","BASE_NUM_DUPL","CONST_RATE","CONST_RATE_DEMI","CONST_RATE_DEMI_EST","CONST_RATE_NO_DUPL"] # all models
counts_file = ".counts_edit"

with open (filename, "r") as genera:
	with open("/groups/itay_mayrose/annarice/model_adequacy/genera/summary_nine_generaX.csv", "w+") as writeFile:
		header = ["Genus", "Model"]
		writer = csv.writer(writeFile)
		writer.writerow(header)

		for genus in genera:
			genus = genus.strip()
			if models_flag == "ALL":
				models = all_models
			if models_flag == "BEST":
				results_sum = "/groups/itay_mayrose/annarice/model_adequacy/genera/" + genus + "/result_sum"
				models = best_model.get_best_model(results_sum)
				models = models.split(",") # returns a list
			if models_flag == "OTHERS":
				models = list(set(all_models) - set(models)) # run on all models that are not the best model
			if models_flag == "DEFINED":
				models = models_defined.split(",")  # returns a list

			for model in models:
				if sanity == 1:
					for i in range(50):  # NEED THIS FOR MA ON SANITY
						name = genus + str(i) + "." + model
						wd = "/groups/itay_mayrose/annarice/model_adequacy/sanity/" + genus + "/CONST_RATE/adequacy_test/" + str(i) + "/"
						co = "/groups/itay_mayrose/annarice/model_adequacy/sanity/" + genus + "/CONST_RATE/adequacy_test/" + str(i) + "/counts.txt"
						tree = "/groups/itay_mayrose/annarice/model_adequacy/sanity/" + genus + "/tree_1"
						cmd = "-c " + wd + " -m " + model + " -id " + genus + " -nt 1 -ns " + str(nsims) + " -ce /groups/itay_mayrose/itaymay/code/chromEvol/chromEvol_source-current/chromEvol -co " + co + " -t " + tree + " -s " + str(sanity) + " -r " + str(results_flag)
						os.system("python ~/create_sh.py " + "-name " + name + " -l " + lang + " -s " + script + " -p \'" + cmd + "\'")
				else:
					name = genus + "." + model
					wd = "/groups/itay_mayrose/annarice/model_adequacy/genera/" + genus + "/" + model + "/"
					co = "/groups/itay_mayrose/annarice/model_adequacy/genera/" + genus + "/" + genus + counts_file
					tree = "/groups/itay_mayrose/annarice/model_adequacy/genera/" + genus + "/tree_1"
					cmd = "-c " + wd + " -m " + model + " -id " + genus + " -nt 1 -ns " + str(nsims) + " -ce /bioseq/chromEvol/chromEvol.exe -co " + co + " -t " + tree + " -s " + str(sanity) + " -r " + str(results_flag)
					row = [genus, models]
					writer.writerow(row)
					os.system("python ~/create_sh.py " + "-name " + name + " -l " + lang + " -s " + script + " -p \'" + cmd + "\'")



# in case of a list of genera --> lst = ["A","B","C"]
'''
if sanity=="1":
	for genus in lst:
		for model in models:
			for i in range(1): # NEED THIS FOR MA ON SANITY
				name = genus + str(i) + "." + model
				wd = "/groups/itay_mayrose/annarice/model_adequacy/sanity/" + genus + "/CONST_RATE/adequacy_test/" + str(i) + "/"
				co = "/groups/itay_mayrose/annarice/model_adequacy/sanity/" + genus + "/CONST_RATE/adequacy_test/" + str(i) + "/counts.txt"
				tree = "/groups/itay_mayrose/annarice/model_adequacy/sanity/" + genus + "/tree_1"
				cmd = "-c " + wd + " -m " + model + " -id " + genus + " -nt 1 -ns " + str(nsims) + " -ce /bioseq/chromEvol/chromEvol.exe -co " + co + " -t " + tree
				os.system("python ~/create_sh.py " + "-name " + name + " -l " + lang + " -s " + script + " -p \'" + cmd + "\'")
else:
	for genus in lst:
		for model in models:
			name = genus + "." + model
			wd = "/groups/itay_mayrose/annarice/model_adequacy/genera/" + genus + "/" + model + "/"
			co = "/groups/itay_mayrose/annarice/model_adequacy/genera/" + genus + "/" + genus + counts_file
			tree = "/groups/itay_mayrose/annarice/model_adequacy/genera/" + genus + "/tree_1"
			cmd = "-c " + wd + " -m " + model + " -id " + genus + " -nt 1 -ns " + str(nsims) + " -ce /bioseq/chromEvol/chromEvol.exe -co " + co + " -t " + tree
			os.system("python ~/create_sh.py " + "-name " + name + " -l " + lang + " -s " + script + " -p \'" + cmd + "\'")
'''