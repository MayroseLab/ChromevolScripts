"""
	Collect all Medians from CCDB per genus and choose the value according to:
	merge counts of duplicate accepted,
	merge counts of duplicate unresolved,
	filter out unresolved when accepted exists
"""

__author__ = 'MichalDrori'

DEBUG_FLAG = 0

import argparse
import csv
from collections import Counter
import sqlite3
from Bio import Phylo
import unicodedata
import pickle
import sys
import os
from collections import namedtuple

CCDB_DATABASE = '/bioseq/chromEvol/CCDB_1.58.db'
#CCDB_DATABASE = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/CCDB_1.46.db'

#--------------------------------------------------def-------------------------------------------------------
def turn_list_to_DB_list(list_name):

	db_names_str=''
	for item in list_name:
		name_lowercase = str(item).lower()
		item_add="'"+name_lowercase+"',"
		db_names_str+=item_add
	final_str = db_names_str[:-1]
	return final_str
#--------------------------------------------------def-------------------------------------------------------
def get_leaf_names(NewickTree):
	spc_names_on_Tree=[]
	tree = Phylo.read(NewickTree, "newick")
	for leaf in tree.get_terminals():
		spc_names_on_Tree.append(leaf.name)
	return spc_names_on_Tree
#--------------------------------------------------def-------------------------------------------------------
def new_90prec_max(last_dist_str):
#example 24=0.270000_16=0.440000_8=0.100000_20=0.120000_28=0.010000_40=0.010000_12=0.030000_32=0.010000_23=0.010000

	if DEBUG_FLAG == 1: print("new_90prec_max------------------------------")
	count_prec_dict=dict()
	split_data = last_dist_str.strip().split('_')
	if DEBUG_FLAG == 1: print('pppppppppppppppppppppppppp')
	if DEBUG_FLAG == 1: print(split_data)
	for item in split_data:
		if DEBUG_FLAG == 1: print('ppppppppppppppppppppppp')
		if DEBUG_FLAG == 1: print(item)
		count = int(item.split('=')[0])
		prec = item.split('=')[1]
		count_prec_dict[count]=prec


	count_max_sort = sorted(count_prec_dict.keys(),reverse=True)
	if DEBUG_FLAG == 1: print(count_max_sort)
	sum_prec=0.0000
	for count in count_max_sort:
		sum_prec+=float(count_prec_dict[count])
		if DEBUG_FLAG == 1: print(sum_prec)
		if DEBUG_FLAG == 1: print(count)
		if sum_prec > 0.1:
			return count
		else:
			print("Still....checking")

	return New_max
#--------------------------------------------------def-------------------------------------------------------
def break_doubles(species_dist_dict):

	clean_list=[]
	for item in species_dist_dict:
		if not item:
			continue
		elif ',' in item:
			# Split counts in case more than one:
			for count_item in item.split(','):
				clean_list.append(count_item.strip())
		else:
			clean_list.append(item)

	return clean_list
#--------------------------------------------------def-------------------------------------------------------
def check_max_value(dist_dict):
	if DEBUG_FLAG == 1: print(dist_dict.keys())
	max_count = max(dist_dict.keys())
	if dist_dict[max_count] >= 0.1:
		return max_count
	else:
		del dist_dict[max_count]
		new_max_count = max(dist_dict.keys())
		return new_max_count
# --------------------------------------------------def-------------------------------------------------------
def calc_distribution(species_dist_dict,path):

	#split data which include list of counts for one count
	min_f = open(path + '_min', 'w')
	max_f = open(path + '_max', 'w')
	dist_f = open(path + '_dist', 'w')

	lines_max=''
	lines_min=''
	lines_dist=''
	for specie in species_dist_dict.keys():
		specie_ = escape_organism_name(specie)
		lines_dist += '>%s\n'%specie_.replace(' ','_')
		numberOfcounts = 0
		Sum_counts = 0
		max_count = -1
		min_count = 999999
		list_single_data = break_doubles(species_dist_dict[specie])
		dist_counts = Counter(list_single_data)
		for count in list_single_data:
			if not count:
				continue
			else:
				Sum_counts+=int(count)
				if int(count) > max_count:
					max_count = int(count)
				if int(count) < min_count:
					min_count = int(count)
				numberOfcounts+=1
		dist_str=''
		dist_dict = dict()
		max_prec=-1
		max_item='null'
		sum_prec=0
		for item in dist_counts.keys():
			dist_str+=str(item)
			prectg = round(float(dist_counts[item])/float(numberOfcounts),2)
			sum_prec+=prectg
			if prectg > max_prec:
				max_prec = prectg
				max_item = item
			dist_str += "=%f_" % prectg
			dist_dict[int(item)]=prectg
		if float(sum_prec) == float(1):
			print('Sum is 1')
		else:
			diff_prec = float(1)-sum_prec
			if DEBUG_FLAG == 1: print("diff_prec %f" %diff_prec)
			if DEBUG_FLAG == 1: print("max_prec Before %f" % max_prec)
			max_prec+=diff_prec
			if DEBUG_FLAG == 1: print("max_prec After %f" % max_prec)
			dist_str = ''
			for item in dist_counts.keys():
				dist_str+=str(item)
				prectg = round(float(dist_counts[item])/float(numberOfcounts),2)
				if item == max_item:
					prectg=max_prec
				dist_str += "=%f_" % prectg
				dist_dict[int(item)] = prectg
		#check max value to be at least 0.1 of dist, otherwise replace with next value:
		#if dist_dict.keys():
		#	final_max = check_max_value(dist_dict)
		#else:
		#	final_max = max_count
		last_dist_str = dist_str[:-1]
		if last_dist_str:
			New_max = new_90prec_max(last_dist_str)
		else:
			New_max=max_count
		#lines_max += '%s\n'%str(final_max)
		if New_max == -1:
			print("%s, No maximum found" %specie_)
		else:
			lines_max += '>%s\n' % specie_.replace(' ', '_')
			lines_max += '%s\n'%str(New_max)
		if min_count == 999999:
			print("%s, No minimum found" % specie_)
		else:
			lines_min += '>%s\n' % specie_.replace(' ', '_')
			lines_min += '%s\n'%str(min_count)
		lines_dist += '%s\n'%str(last_dist_str)

	max_f.write(lines_max)
	min_f.write(lines_min)
	dist_f.write(lines_dist)

	return
#----------------------------------------------------------------------------------
def escape_organism_name(organism_name):
	if organism_name is None:
		return None


	escaped_organism = organism_name
	#Added space after ï¿½ to be aligned with the name in the plant list:
	escaped_organism = escaped_organism.replace("Ã—", "x")
	escaped_organism = escaped_organism.replace(" ", "_")
	escaped_organism = escaped_organism.replace(",", "_")
	escaped_organism = escaped_organism.replace("-", "_")
	escaped_organism = escaped_organism.replace("'", "_")
	escaped_organism = escaped_organism.replace("/", "_")
	escaped_organism = escaped_organism.replace("&", "AND")
	escaped_organism = unicodedata.normalize('NFKD', escaped_organism).encode('ascii','ignore')
	escaped_organism = str(escaped_organism, encoding='ascii')

	return escaped_organism
#--------------------------------------------------------------------------------------
def get_species_list_for_genus (genus,db_curser):
	if DEBUG_FLAG == 1: print("##" + genus + "\n")
	if DEBUG_FLAG == 1: print("-----------------------------\n")
	species_names=[]
	species_list_unique=[]
	query_genus_species_db = "SELECT resolved_name_full FROM resolved_chrom_counts WHERE genus like '%s'" % genus
	db_curser.execute(query_genus_species_db)
	rows_genus_db = db_curser.fetchall()
	for line in rows_genus_db:
		species_names.append(line[0])
	species_list_unique=set(species_names)
	if DEBUG_FLAG == 1: print(species_list_unique)
	return species_list_unique
#--------------------------------------------------------------------------------------
def get_counts_for_species (genus,db_curser):

	query_species_counts_db = "SELECT id,genus,resolved_name_full,resolved_name,median,status FROM resolved_chrom_counts WHERE resolved_name like '%s%%'" % genus
	db_curser.execute(query_species_counts_db)
	rows_genus_db = db_curser.fetchall()

	return rows_genus_db
#--------------------------------------------------------------------------------------
def get_ploidy_data_all(ploidy_file,file_handle):
	ploidy_dict = {}
	chromCount_dict = {}
	with open(ploidy_file,'r') as f:
		reader = csv.DictReader(f)
		indx=0
		for taxon_line in reader:
			species_name=taxon_line['Taxon']
			if DEBUG_FLAG == 1: print(">>>>>>>>>>>>>>>>>>>>>>>>>>>")
			if DEBUG_FLAG == 1: print(species_name)
			chromCount_dict[species_name] = taxon_line['Chromosome count']
			ploidy_dict[species_name] = taxon_line['Ploidy inference']
			indx+=1
			if DEBUG_FLAG == 1:
				file_handle.write("%s,%s,%s\n" % (species_name,chromCount_dict[species_name],ploidy_dict[species_name]))
		return (ploidy_dict,chromCount_dict)
#--------------------------------------------------------------------------------------
def recalc_median(species_name,db_curser):
	list_parsed_n=[]
	query_species_parsed_n_db = "SELECT parsed_n FROM resolved_chrom_counts WHERE resolved_name LIKE '%s'" % species_name
	db_curser.execute(query_species_parsed_n_db)
	rows_genus_db = db_curser.fetchall()
	for line in rows_genus_db:
		line_list=list(line)
		if "None" in str(line_list):
			continue
		else:
			for item in line_list:
				parsed_n = item.replace(" ","")
				parsed_n = parsed_n.split(",")
				if '' in parsed_n:
					continue
				else:
					for parsed_val in parsed_n:
						list_parsed_n.append(int(parsed_val))
	if not list_parsed_n:
		new_median=0
	else:
		new_median = int(statistics.median(list_parsed_n))
	return new_median
#--------------------------------------------------------------------------------------
def break_comma_db_rows(db_rows): 
	new_db_rows = []
	for line in db_rows:
		if not line[1]:
			continue
		if ',' not in line[1]:
			new_db_rows.append(line)
			continue
		vals = line[1].split(',')
		for v in vals:
			new_line = []; 
			new_line.append(line[0])
			new_line.append(v.strip())
			new_db_rows.append(new_line)
	return new_db_rows
#--------------------------------------------------------------------------------------
def get_db_rows_value(db_rows, flag): 

	prev_spec = ''
	n_list = []
	n_counts = []
	counts = {}
	counts_list = []
	total_counts = 0
	results = []
	
	for line in db_rows:
		cur_spec = line[0]
		
		if cur_spec != prev_spec and prev_spec != '':
		
			# append species group to results
			sorted_l = sorted(n_list)
			counts_list = []
			u_score = ''
			dist_str = ''
			n = 0
			len_list = len(n_list)
			sum = 0
			for i in sorted_l: 
				counts_list.append(counts[str(i)])
				freq = counts[str(i)] / total_counts
				if n < len_list-1:
					ffreq = "{:.2f}".format(freq)
				else:
					ffreq = "{:.2f}".format(round(1-sum,2))
				sum = sum + float(ffreq)
				dist_str += f'{u_score}{i}={ffreq}'
				u_score = '_'
				n += 1
			
			if flag == 'max':
				val = max(n_list)
			elif flag == 'min':
				val = min(n_list)
			elif len(n_list) > 1:
				val = dist_str
			else:
				val = n_list[0]
			
			new_line = []; 
			new_line.append(prev_spec)
			new_line.append(val)
			results.append(new_line)
			
			# reset vars
			prev_spec = cur_spec
			n_list = []
			n_counts = []
			counts = {}
			total_counts = 0
		
		# on first loop
		if prev_spec == '':
			prev_spec = cur_spec
			
		# add to count lists
		ival = int(line[1])
		cval = line[1]
		if ival not in n_list: 
			n_list.append(ival)
			counts[cval] = 0
		counts[cval] += 1
		total_counts += 1
	
	# last value
	if len(n_list) > 0: 
		sorted_l = sorted(n_list)
		counts_list = []
		u_score = ''
		dist_str = ''
		n = 0
		len_list = len(n_list)
		sum = 0
		for i in sorted_l: 
			counts_list.append(counts[str(i)])
			freq = counts[str(i)] / total_counts
			if n < len_list-1:
				ffreq = "{:.2f}".format(freq)
			else:
				ffreq = "{:.2f}".format(round(1-sum,2))
			sum = sum + float(ffreq)
			dist_str += f'{u_score}{i}={ffreq}'
			u_score = '_'
			n += 1
			
		if flag == 'max':
			val = max(n_list)
		elif flag == 'min':
			val = min(n_list)
		elif len(n_list) > 1:
			val = dist_str
		else:
			val = n_list[0]
		
		new_line = []; 
		new_line.append(line[0])
		new_line.append(val)
		results.append(new_line)
			
		return results

if __name__ == "__main__":

	tree_or_genera = sys.argv[1]
	Output_ForCounts = sys.argv[2] 	# log file
	Tree_or_Name = sys.argv[3] 		# 0 -> Tree file 1 -> list of genera in file
	counts_type = sys.argv[4] 		# median, max, min, distribution
	genera_list=[]
	taxa_list=[]
	out_file = "Counts.txt"
	
	file_out = open(Output_ForCounts + "/COUNTS_DATA_LOG.txt", mode="w", encoding='utf8', newline='')
	if Tree_or_Name == '0':
		species_list = get_leaf_names(tree_or_genera)
		file_out.write("Taxa names on tree:\n")
		for spc in species_list:
			taxa_name = spc.strip().replace('_',' ')
			taxa_list.append(taxa_name)
			file_out.write('%s\n' %taxa_name)
			genus = spc.split('_')[0]
			if genus not in genera_list:
				genera_list.append(genus)
	elif Tree_or_Name == '1':
		with open(tree_or_genera,'r') as f:
			genera_list = [g.strip() for g in f.readlines()]


	species_dist_dict = dict()
	counts_file = Output_ForCounts + '/countsFile'
	counts_file = open(counts_file, 'w')
	#species_names = []
	#species_median_dict = {}
	#species_status_dict = {}
	#species_parsed_n_dict = {}
	
	# db connection
	if not os.path.exists(CCDB_DATABASE): 
		file_out.write("file %s does not exist\n" %CCDB_DATABASE)
	conn_db = sqlite3.connect(CCDB_DATABASE)
	db_curser = conn_db.cursor()
	taxa_list_str = turn_list_to_DB_list(taxa_list)
	genera_list_str = turn_list_to_DB_list(genera_list)
	
	# create query string
	if Tree_or_Name != '1':
		if counts_type == 'median': 
			query_taxa_counts_db = "SELECT distinct resolved_name,median from resolved_chrom_counts where lower(resolved_name) in (%s)" % taxa_list_str
		else: 
			query_taxa_counts_db = "SELECT resolved_name, parsed_n from resolved_chrom_counts where lower(resolved_name) in (%s) ORDER by resolved_name" % taxa_list_str
	else:
		query_taxa_counts_db = "SELECT distinct resolved_name,median from resolved_chrom_counts where lower(genus) in (%s)" % genera_list_str
	
	# execute query string
	file_out.write('%s\n' % query_taxa_counts_db)
	db_curser.execute(query_taxa_counts_db)
	rows_genus_db = db_curser.fetchall()
	
	# derive counts according to counts type
	if counts_type == 'median': 
		results = rows_genus_db
	else: # min, max, distribution
		new_rows_genus_db = break_comma_db_rows (rows_genus_db)
		file_out.write("breaking down rows:\n")
		for line in new_rows_genus_db:
			file_out.write("{line}\n")
		results = get_db_rows_value(new_rows_genus_db, counts_type)
		
	# write to counts file
	count_counts = 0
	for line in results:
		specie_name = line[0]
		if line[1]:
			val = line[1]
		else:
			val = 'x'
		counts_file.write('>%s\n' % specie_name.replace(' ','_'))
		counts_file.write('%s\n' % val)
		count_counts += 1
	
	# write species list for OneTwoTree
	species_list = Output_ForCounts + '/species_list'
	f_species = open(species_list, 'w')
	comma = ""
	for line in results:
		specie_name = line[0]
		f_species.write (comma+specie_name)
		comma = ','
	f_species.close()
		
	# Check number of counts:
	if count_counts < 5:
		less_f = open(Output_ForCounts + '/LESS_THAN_5_COUNTS', 'w')
		less_f.close()
	counts_file.close()

