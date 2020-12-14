from MA_defs import *
from MA_utils import *
from data_processing import process_data
from analysis import get_stats


def create_simulated_stats_distribution(out_dir,nsims,main_res_dir): #[v, e, r, u]
    simulated_counts_stats_dist = []
    for i in range(nsims):  # for each simulation
        sim_counts = process_data.get_counts(out_dir + str(i) + "/simCounts.txt")
        print(main_res_dir)
        simulated_counts_statistics = get_stats.calculate_statistics(sim_counts,out_dir + str(i) + "/simStats", main_res_dir + "/mlAncestors.tree",out_dir + str(i) + "/simCounts.txt")
        simulated_counts_stats_dist.append(simulated_counts_statistics)
    return simulated_counts_stats_dist


def handle_distributions(dist,stats,res,filename1,filename2):
    # print distribution to file
    with open(filename1, "w") as distribution_file:
        for i in range (len(stats)):
            sim_stat_dist = [round(x[i],4) for x in dist]
            distribution_file.write(str(sim_stat_dist))
            distribution_file.write("\n")
    with open(filename2,"w") as dist_vec:
        dist_vec.write(str(res))

def test_adequacy(sim_stats, stats,filename):
    '''
    go over each statistic in the list and see if is adequate
    '''
    adequacy_lst = []
    stat_lst = []
    with open(filename, "w+") as percentiles:
        for i in range(len(stats)):
            sim_stat_dist = [x[i] for x in sim_stats]  # simulated ith distribution, x the item in each simulated list
            stat_star_upper = np.percentile(sim_stat_dist, 97.5)  # calculate the upper limit
            stat_star_lower = np.percentile(sim_stat_dist, 2.5)  # calculate the lower limit
            model = 0
            if stats[i] <= stat_star_upper and stats[i] >= stat_star_lower:
                model = 1
            adequacy_lst.append(model)
            stat_lst.append(stats[i])
            percentiles.write(str(round(stat_star_lower, 4)) + "," + str(round(stat_star_upper, 4)) + "\n")
    return (adequacy_lst,stat_lst)

def post_analysis(adequacy_results,stats_results,model_name,filename,id):
    CE_res_filename, expectation_file, mlAncTree, root_freq_filename, sim_control, statistics_names = fixed_vars()
    with open(filename, "w+") as results_file:
        if all(x==1 for x in adequacy_results):
                print ("In " + id + ", " + model_name + " is adequate for all statistics", file = results_file)
                for i in range(len(stats_results)):
                    print (str(statistics_names[i] + " = " + str(stats_results[i])), file = results_file)
                return adequacy_results
        print("In " + id + ", " + model_name + " is: ", file = results_file)
        for i in range(len(adequacy_results)):
            if adequacy_results[i]==0:
                    print ("Not adequate for " + str(statistics_names[i]), file = results_file)
            if adequacy_results[i]==1:
                    print ("Adequate for " + str(statistics_names[i]), file = results_file)
        #for i in range(len(stats_results)):
         #   print (str(statistics_names[i] + " = " + str(stats_results[i])), file = results_file)
        return adequacy_results

def model_adequacy(out_dir,orig_counts_stats, model_name,max_for_sims,nsims,id,main_res_dir):
    sim_dist = create_simulated_stats_distribution(out_dir,nsims,main_res_dir)
    adequacy_results,stats_results = test_adequacy(sim_dist, orig_counts_stats,out_dir + "percentiles_" + model_name)
    handle_distributions(sim_dist, orig_counts_stats,adequacy_results, out_dir + "stats_dist_sims", out_dir + "adequacy_vec")
    results_lst = post_analysis(adequacy_results,stats_results, model_name, out_dir + model_name + "_MA.res",id)
    l = []
    for i in range(nsims):
            l.append(str(i))
    targz_dir(out_dir, l, "zipped.tar.gz", True)
    return results_lst