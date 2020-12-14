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
#from taxonome.taxa.file_csv import load_taxa
#from taxonome.taxa.name_selector import NameSelector
#from taxonome.tracker import CSVTaxaTracker, _flatten_name, CSVListMatches
import unicodedata


import pickle
import sys
import os
from collections import namedtuple

# Running line:
#python /groups/itay_mayrose/michaldrori/Pipeline_Files/CCDB/make_ccdb_counts_files_for_chromevol.py Pipeline_Files/Alignment_TaxonList/Taxon_list_origin.txt

CCDB_DATABASE = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/CCDB_1.46.db'


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


def check_max_value(dist_dict):
	if DEBUG_FLAG == 1: print(dist_dict.keys())
	max_count = max(dist_dict.keys())
	if dist_dict[max_count] >= 0.1:
		return max_count
	else:
		del dist_dict[max_count]
		new_max_count = max(dist_dict.keys())
		return new_max_count


#max,dist = calc_distribution(species_dist_dict)
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
	#Added space after � to be aligned with the name in the plant list:
	escaped_organism = escaped_organism.replace("×", "x")
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

#D:\Michal\scripts\make_ccdb_counts_files_for_chromevol.py

if __name__ == "__main__":


	#Count_type = 'Min' #options are: median, max, dist (distribution)
	genera_list_var = sys.argv[1]
	Output_ForCounts = sys.argv[2] #'/groups/itay_mayrose/michaldrori/MD_ChromEvol/Counts_ccdb'
	File_or_Name = sys.argv[3] # 0 -> File with a list of genera or 1 -> list of genera in string 'AAA,BBB,CCC'
	#species_list=[]
	out_file = "Counts.txt"
	#out_handle = open(out_file,mode="w",encoding='utf8',newline='')
	#writer = csv.writer(out_handle,delimiter=',')
	#writer.writerow(['genus','resolved_name_full','resolved_name','median','status'])
	#with open('/groups/itay_mayrose/michaldrori/Pipeline_Files/Alignment_TaxonList/Angio_Missing.txt','r') as f:
	if File_or_Name == '0':
		with open(genera_list_var,'r') as f:
			genera_list = [g.strip() for g in f.readlines()]
	elif File_or_Name == '1':
		genera_list= genera_list_var.split(',')

	file_out=open(Output_ForCounts + "/COUNTS_DATA_LOG.txt",mode="w",encoding='utf8',newline='')
	file_LessThan5Counts=open(Output_ForCounts + "/Genera_withLessThan5Counts.txt",mode="w",encoding='utf8',newline='')
	for genus in genera_list:
		species_dist_dict=dict()
		count_counts=0
		print(genus)
		counts_file = Output_ForCounts + '/' + genus + '.counts'
		counts_file=open(counts_file,'w')
		species_names=[]
		species_median_dict={}
		species_status_dict={}
		species_parsed_n_dict={}
		conn_db = sqlite3.connect(CCDB_DATABASE)
		db_curser = conn_db.cursor()
		query_genus_species_db = "SELECT id,genus,resolved_name_full,resolved_name,parsed_n,median,status from resolved_chrom_counts WHERE genus is '%s'" % genus
		#query_genus_species_db = "SELECT resolved_name FROM resolved_chrom_counts WHERE genus is '%s'" % genus
		db_curser.execute(query_genus_species_db)
		rows_genus_db = db_curser.fetchall()
		for line in rows_genus_db:
			specie_name = line[3]
			parsed_n = line[4]
			median_val = line[5]
			status = line[6]
			#Check for duplicate species names:
			if specie_name in species_names:
				if DEBUG_FLAG == 1: print(specie_name)
				if DEBUG_FLAG == 1: print(parsed_n)
				if DEBUG_FLAG == 1: print(status)
				#create data arrayes for max/dist calc
				#if status: # == 'Accepted':
				if specie_name in species_dist_dict.keys():
					species_dist_dict[specie_name].append(parsed_n)
				else:
					species_dist_dict[specie_name]=[parsed_n]
				#Check if same median value:
				if median_val == species_median_dict[specie_name]:
					#No need to add data
					continue
				#Check if different median and different status:
				else:
					#Check if Different status and take Accepted median:
					if status != species_median_dict[specie_name]:
						#Update values according to Accepted specie:
						if status == 'Accepted':
							#First remove the existing value if exists:
							if specie_name in species_names:
								file_out.write("%s is being removed, Accepted median was found\n" %specie_name)
								species_names.remove(specie_name)
								del species_median_dict[specie_name]
								del species_status_dict[specie_name]
								del species_parsed_n_dict[specie_name]
							species_names.append(specie_name)
							species_median_dict[specie_name]=median_val
							species_status_dict[specie_name]=status
							species_parsed_n_dict[specie_name]=parsed_n
							if specie_name in species_dist_dict.keys():
								species_dist_dict[specie_name].append(parsed_n)
							else:
								species_dist_dict[specie_name] = [parsed_n]
						else:
							continue
					else:
						#Different values to same status:
						file_out.write('## Need to re-calc median for %s\n' % specie_name)
						continue
			#Specie name not in list:
			else:
				species_names.append(specie_name)
				species_median_dict[specie_name]=median_val
				species_status_dict[specie_name]=status
				species_parsed_n_dict[specie_name]=parsed_n
				if specie_name in species_dist_dict.keys():
					species_dist_dict[specie_name].append(parsed_n)
				else:
					species_dist_dict[specie_name]=[parsed_n]

		if DEBUG_FLAG == 1: print("DDDDDDDDDDDDDD - > Dist by name:")
		if DEBUG_FLAG == 1: print(species_dist_dict)
		#Calc dist and max for each
		calc_distribution(species_dist_dict,Output_ForCounts + '/' + genus + '.counts')

		#species_list_unique=set(species_names)
		#file_out.write(str(species_list_unique))
		#file_out.write(str(len(species_list_unique)))
		for specie in species_names:
			if DEBUG_FLAG == 1: print(specie)
			if DEBUG_FLAG == 1: print("%s,%s\n" %(genus,specie))
			if DEBUG_FLAG == 1: file_out.write("%s,%s,%s,%s,%s\n" % (genus, specie,species_median_dict[specie],species_parsed_n_dict[specie],species_status_dict[specie]))
			specie_=escape_organism_name(specie)
			#specie_=specie_x.replace(' ','_')
			if not species_median_dict[specie]:
				file_out.write('### None value for  %s' % specie_name)
				continue
			else:
				counts_file.write('>%s\n' % specie_)
				counts_file.write('%s\n' % species_median_dict[specie])
				count_counts+=1
		#Check number of counts:
		if count_counts < 5:
			file_LessThan5Counts.write(genus+'\n')

