from MA_defs import *
import os
from ete3 import Tree
import re
import argparse

def calculate_statistics(counts,filename, tree_file,simulated_counts_file = False):
    '''
        given list of counts produces statistics: variance,min,max,entropy
        ########## ADD MP OF NUMBER OF TRANSITIONS
    :param counts: list of chromosome numbers
    :return: list of statistics representing the counts
    '''
    # variance
    v = np.var(counts)

    # range = max - min
    r = max(counts) - min(counts)

    # enthropy, calculates the probabilities
    d = {}
    for i in counts:
        d[i] = counts.count(i)
    prob_lst = [x / len(counts) for x in list(d.values())]
    e = sc.entropy(prob_lst)

    # unique counts
    counts_set = set(counts)
    u = len(counts_set)

    # parsimony
    p = fitch(tree_file,simulated_counts_file)

    with open(filename, "w+") as stats:
        stats.write(str(round(v,4)) + "," + str(round(e,4)) + "," + str(round(r,4)) + "," + str(round(u,4)) + "," + str(round(p,4)))

    return ([v, e, r, u, p])

def fitch (tree_file, c = False):
    print(tree_file)
    t = Tree(tree_file, format = 1)
    score = 0

    if c: # if there's a counts file, the analysis is of a simulated dataset
        print(c)
        d = {}
        with open(c, "r") as counts:
            for line in counts:
                line = line.strip()
                if line.startswith(">"):
                    key = line[1:]
                else:
                    val = line
                    d[key] = val

    for node in t.traverse("postorder"):
        if not node.is_leaf():  # internal node
            lst = []  # list version
            intersect, union = None, None
            for child in node.get_children():
                if child.is_leaf():  # if the child is a tip - parse number from tip label
                    if c:  # if there is a dictionary --> take the number from it --> the tree is a simulated tree
                        name = re.search("(.*)\-\d+", child.name)
                        if name:
                            num = {int(d.get(name.group(1)))}
                    else:  # the tree is the original tree
                        tmp = re.search("(\d+)", child.name)
                        if tmp:  # there is a number at the tip, and not X
                            num = {int(tmp.group(1))}
                else:  # if the child is an internal node - take number
                    num = child.name
                lst.append(num)
            intersect = lst[0] & lst[1]
            union = lst[0] | lst[1]
            #print("intersect " + str(intersect))
            #print("union " + str(union))
            if len(intersect) == 0:
                result = union
                score += 1
            else:
                result = intersect
            node.name = result

    return(score)