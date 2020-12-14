#from defs import *

def get_counts(filename):
    '''
        reads the .counts_edit file and extracts the counts
    :param filename: supplied by the user
    :return: list of counts
    '''
    with open(filename, "r") as tmp_counts_file:
        counts = []
        for line in tmp_counts_file:
            line = line.strip()
            if line.startswith('>'):
                continue
            else:
                if line=="x":
                    continue
                counts.append(int(line))
    return (counts)

'''
def extract_empirical_data(counts_lst):
        get the min chr in data, min chr allowed, max chr in data, and max chr allowed, to be printed to the control file
        ################### CURRENTLY- MIN ALLOWED IS -1 AND MAX ALLOWED IS *10
        empirical data can also be extracted from the CE results file
    :param counts_lst: get_counts output
    :return: list of original counts min and max variations

    return([min(counts_lst),max(counts_lst),min(counts_lst)-1,max(counts_lst)+10])
    
'''
