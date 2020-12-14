import os

from utils import *
from defs import *
from data_processing import process_data
from data_processing import best_model
from data_processing import simulations
from analysis import get_stats
from analysis import test_adequacy

if __name__ == '__main__':
    id, main_res_dir, in_model, num_of_trees, sims_per_tree, CE_path, counts_file, tree_full_path = get_arguments()
    CE_res_filename, expectation_file, mlAncTree, root_freq_filename, sim_control, statistics_names = fixed_vars()
    m = len(in_model)
    for k in range(m):  # run over all models or a single model
        model = in_model[k]
        for i in range(num_of_trees):
            ind = str(i + 1)
            output_dir = main_res_dir + "/adequacy_test/"
            if not os.path.exists(output_dir):
                res = os.system(
                    "mkdir -p " + output_dir)  # -p allows recusive mkdir in case one of the upper directories doesn't exist
            original_counts = process_data.get_counts(counts_file)
            original_counts_statistics = get_stats.calculate_statistics(original_counts, output_dir + "orig_stats",
                                                                        main_res_dir + mlAncTree)
            max_for_simulations = simulations.run_MA(main_res_dir, output_dir + sim_control, main_res_dir,
                                                     output_dir, original_counts, model, sims_per_tree,
                                                     tree_full_path, CE_path)
            adequacy_lst = test_adequacy.model_adequacy(output_dir, original_counts_statistics, model,
                                                        max_for_simulations, sims_per_tree, id, main_res_dir)

        # id - job num/name
        # main_res_dir - chromevol results path for spesicfic model BASE_NUM/
        # in_model - the model the user want to check or selected one
        # num_of_trees - for now only 1 tree
        # sims_per_tree - num of sim (1000)
        # CE_path - chromevol exe path
        # counts_file -
        # tree_full_path -
        #
