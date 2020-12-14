from MA_defs import *
#import pandas as pd

def get_best_model(filename):
    data = pd.read_csv(filename, sep="\t", header=None)
    tmp = data.loc[data[3] == 0,0].values[0]
    return tmp # the name of the best model

#def create_freq_file(line,freqs,upper_bound):
def create_freq_file(line, freqs):
    tmp = line
    # write the tmp to a file, with the aid of emp_data
    with open(freqs, "w+") as root_freq:
        text = tmp.split()
        first = re.search("F\[(\d+)\]", text[0])
        first = int(first.group(1))
        last = re.search("F\[(\d+)\]", text[-1])
        last = int(last.group(1))
        #for i in range(1, first):
        #    print("F[" + str(i) + "]=0", file=root_freq)
        for i in range(0, last - first + 1):
            print(text[i], file=root_freq)
        #for i in range(last + 1, upper_bound):
         #   print("F[" + str(i) + "]=0", file=root_freq)

#def get_params(filename, freqs, upper_bound):
def get_params(filename, freqs):
    '''
        produce a dictionary of parameters
    '''
    params_dict = {}
    line_cntr = 0
    with open(filename, "r") as params_file:
        for line in params_file:
            line_cntr = line_cntr + 1
            if line_cntr == 6:  # reached the parameters part
                line = line.strip()
                tree_length = re.search("(\d+)", line).group(1)
                params_dict["_simulationsTreeLength"] = tree_length
            if line_cntr > 15: # reached the parameters part
                line = line.strip()
                tmp = re.search("(^[^F].*)\s(.*)",line) # the line doesn't start with F, indicating the root frequencies results
                if tmp:
                    key = tmp.group(1)
                    val = int(tmp.group(2)) if key=="DysBnumBER" else float(tmp.group(2))
                    params_dict[key] = val # key = name of parameter, val = parameter's value
                else: # reached the root frequencies part
                    #create_freq_file(line,freqs,upper_bound)
                    create_freq_file(line, freqs)
                    break

    d = {}
    d = {"LOSS_CONST": "_lossConstR", "GAIN_CONST": "_gainConstR", "DUPL": "_duplConstR","DysBnumBER_R": "_baseNumberR", "DysBnumBER": "_baseNumber", "HALF_DUPL":"_demiPloidyR"}
    for key in d:
        if key in params_dict:
            params_dict[key]
            params_dict[d[key]] = params_dict[key]
            del params_dict[key]
    return (params_dict)
