__author__ = 'Lior'

"""
	Script:

	Description: Uses ploiDB and sexual system data from CSV file to create a
					traits data file for BayesTraits, following the format required
					by the program.

	Inputs: ploiDB, sexual system CSV

	Outputs: traits data file

	Parameters:
"""

import csv
import argparse

def get_sexual_system_data(ss_file, cat_type):
	"""
	Creates a dictionary of the sexual system data
	cat_type - sexual system categorization type:
		1 - Dio vs Herm, xDio, xMono
		2 - Dio & xDio vs herm & xMono
	"""
	assert cat_type == 1 or cat_type == 2, "Unknown sexual categorization type (either 1 or 2)"
	ss_dict = {}
	# define sexual system categories
	if cat_type == 1:
		category0 = ['Herm','xMono','xDio']
		category1 = ['Dio']
	elif cat_type == 2:
		category0 = ['Herm','xMono']
		category1 = ['Dio','xDio']
	# parse csv
	with open(ss_file,'r') as mfh:
		csvReader = csv.reader(mfh)
		next(csvReader)
		for row in csvReader:
			if row[1] == 'yes':	# only take species which appear in trees
				species = row[0]
				ss_raw = row[3] # may be empty or contain multiple options, separated by '|'
				if ss_raw == '':
					ss_dict[species] = '-'
				else:
					ss_list = ss_raw.split('|')
					categories = set()	# which categories (0/1) were found for species
					for ss in ss_list:
						if ss in category0:
							categories.add(0)
						elif ss in category1:
							categories.add(1)
						else:
							print("unknown sexual system",ss)
							raise RuntimeError	# unknown sexual system
					# if all sexual systems belong to the same category, take this category
					if len(categories) == 1:
						ss_dict[species] = str(categories.pop())
					elif len(categories) == 2:
						ss_dict[species] = '-'
					else:
						print("0 sexual system categories or more than 2 were found")
						raise RuntimeError	# 0 categories or more than 2 were found
	return ss_dict

def get_ploidy_data(ploidy_file, conf_threshold = 0.95):
	"""
	Creates a dictionary of the ploidy levels of species for a given genus.
	For species without chromosome counts, only those with combined score >= conf_threshold will be taken.
	For species with counts, the same threshold is used, but only on phylogeny robustness score.
	Species below the threshold will be treated as NA data.
	"""
	if not conf_threshold:
		conf_threshold = 0.95
	ploidy_dict = {}
	chromCount_dict = {}
	with open(ploidy_file,'r') as f:
		reader = csv.DictReader(f)
		for taxon_line in reader:
			chromCount_dict[taxon_line['Taxon']] = taxon_line['Chromosome count']
			if taxon_line['Chromosome count'] == 'x':
				if float(taxon_line['Combined score']) >= conf_threshold:
					ploidy_dict[taxon_line['Taxon']] = taxon_line['Ploidy inference']
				else:
					ploidy_dict[taxon_line['Taxon']] = '-'
			else:
				if float(taxon_line['Phylogeny robustness score']) >= conf_threshold:
					ploidy_dict[taxon_line['Taxon']] = taxon_line['Ploidy inference']
				else:
					ploidy_dict[taxon_line['Taxon']] = '-'
		return (ploidy_dict,chromCount_dict)


def combine_data(ss_dict,ploidy_dict):
	"""
	receives the sexual system data dictionary and ploidy dictionary
	and combines them into one dictionary. Raises a warning in case
	the species do not match.
	"""
	ss_species = set(ss_dict.keys())
	ploidy_species = set(ploidy_dict.keys())
	# species in ss dict but not in ploidy dict
	missing_from_ploidy = ss_species - ploidy_species
	if len(missing_from_ploidy) > 0:
		print('Species missing from ploidy data:\n' + list(missing_from_ploidy).join(", "))
	# species in ploidy dict but not in ss dict
	missing_from_ss = ploidy_species - ss_species
	if len(missing_from_ss) > 0:
		print('Species missing from sexual system data:\n' + list(missing_from_ss).join(", "))

	# combine dictionaries
	traits_dict = {}
	for species in ss_dict:
		ss = ss_dict[species]
		if species in missing_from_ploidy:
			ploidy = "-"
		else:
			ploidy = ploidy_dict[species]
		traits_dict[species] = [ss,ploidy]
	# add species with missing ss data
	for species in missing_from_ss:
		traits_dict[species] = ["-",ploidy_dict[species]]

	return traits_dict

def check_variability(traits_dict):
	"""
	Checks if both ploidy and sexual system have at least one 1 and one 0
	"""
	sexual_systems = [x[0] for x in traits_dict.values()]
	ploidy_infers = [x[1] for x in traits_dict.values()]
	if '1' in sexual_systems and '0' in sexual_systems and '1' in ploidy_infers and '0' in ploidy_infers:
		return True
	else:
		return False


def print_traits(traits_dict,out_file):
	"""
	Print the traits in BayesTraits format:
	species   ploidy   sexual system
	"""
	with open(out_file,'w') as fo:
		for species in traits_dict:
			fo.write('\t'.join([species,traits_dict[species][1],traits_dict[species][0]])) #)print('\t'.join([species,traits_dict[species][1],traits_dict[species][0]]),file=fo)

def prepare_traits_file(ploidy_file,ss,output,categorization,ploidy_threshold=None):
	"""
	main function, runs the whole process
	"""
	ss_dict = get_sexual_system_data(ss,categorization)
	ploidy_dict = get_ploidy_data(ploidy_file,ploidy_threshold)
	traits_dict = combine_data(ss_dict, ploidy_dict)
	if check_variability(traits_dict):
		print_traits(traits_dict,output)
		return True
	else:
		return False

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--ploidy_file','-p', help='path to ploidy csv file', required=True)
	parser.add_argument('--ss','-s', help='path to sexual systems csv', required=True)
	parser.add_argument('-threshold','-t',help='Confidence threshold for ploidy inference')
	parser.add_argument('--output','-o', help='output traits file', required=True)
	parser.add_argument('--categorization','-c',type=int,choices=[1,2],help='Sexual system categorization type', required=True)
	args = parser.parse_args()

	prepare_traits_file(args.ploidy_file,args.ss,args.output,args.categorization,args.threshold)