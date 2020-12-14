import sys
import os
from Bio import Phylo
import urllib.request
from urllib.error import HTTPError


def retrieve_genus_counts(genus):
	print(genus)
	req = urllib.request.Request(url=r"http://ccdb.tau.ac.il/services/countsFull/?genus=" + genus + "&format=json")
	try:
		handler = urllib.request.urlopen(req)
	except HTTPError as e:
		content = e.read()
		return (content)


def get_leaf_names(NewickTree):
	spc_names_on_Tree=[]
	tree = Phylo.read(NewickTree, "newick")
	for leaf in tree.get_terminals():
		spc_names_on_Tree.append(leaf.name)
	return spc_names_on_Tree


#NewickTree = sys.argv[1]
NewickTree = '/groups/itay_mayrose/michaldrori/MD_ChromEvol/SERVER_chromEvol/TreeFile_single'

spc_name = get_leaf_names(NewickTree)
print(spc_name)
genus = spc_name[0].split('_')[0]

page = retrieve_genus_counts(genus)
f = open('/groups/itay_mayrose/michaldrori/MD_ChromEvol/SERVER_chromEvol/counts_ccdb.txt', "wb")
f.write(page)
f.close()

