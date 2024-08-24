"""
Explore how the genetic architecture of local adaptation change through time.
Calculate basic statistics for each time sample.
The multiple output tables can be further analyzed in timeSeries_multiTimes.py

# Usage:
python3 timeSeries_slimHistory.py \
     -i /home/tianlin/Documents/github/data/slim_data/glacial_history/historical_optimum0_timeSeries_4PolyLevels/M2b_smallHighVm_lowMig_clineMap/Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-07_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed342902065499469520_tick

tianlin.duan42@gmail.com
2024.06.05
"""
############################# modules #########################################
import msprime
import tskit
import pyslim
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.spatial import distance_matrix
import numpy as np
import random
from time import time
import sys # for sys.exit()
import allel # for allel.weir_cockerham_fst()
import os # mkdir

# custom package
# import handletrees

# ############################# options #############################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
                    help='Input .trees file suffix with its absolute path',
                    type=str)
# parser.add_argument('-p', '--plot',
#                     help='1 for generating plots and 0 for not plotting',
#                     type=int, default=1)
args = parser.parse_args()
############################# functions #######################################
def LF_fitness(ind_x, ind_y, phenotype,
               optima, dist_mate, sigma_w):
    """Takes x and y coordinates (array-like), phenotypes (array-like),
        environmental optima (array-like), a maximum distance for mating (single value),
        and a standard deviation of fitness function (single value).
    Return two arrays of the average fitness of local and foreign individuals, respectively.
    v2: time efficient but not memory efficient.
    Try the other version when the number of individuals is large.
    """
    optima=np.array(optima)
    coord = list(zip(ind_x, ind_y))
    dist_matrix = distance_matrix(coord, coord)
    isLocal_matrix = dist_matrix <= dist_mate
    # Calculate fitness matrix (w_matrix[i][j]: individual j at location i) with broadcast
    w_matrix = 1.0 + stats.norm.pdf(np.array(phenotype),
                                    np.array(optima)[:, np.newaxis],
                                    sigma_w)
    # relative fitness
    w_matrix = w_matrix / np.mean(w_matrix, 1)[:, np.newaxis]
    #Average relative fitness of local and foreign individuals
    w_local = np.mean(np.ma.masked_array(w_matrix, np.invert(isLocal_matrix)), 1)
    w_foreign = np.mean(np.ma.masked_array(w_matrix, isLocal_matrix), 1)
    return w_local, w_foreign

def calculate_lfmut(ts, dist_mate, sigma_w):
    '''
    Calculate the contribution of each mutation to LF
    input: tree sequences, maximum mating distance and the sd of Gaussian
            stabilizing selection;
    return: mean LF and an array of LF_mut
    Version 2: Randomize the distribution of the mutation and skip calculation
    for neutral mutations
    '''
    # Make a list of lists of mutation effects,
    # corresponding to phenotypic effect of allele 0/1/2... at each site
    ind_x, ind_y, ind_z = zip(*ts.individuals_location)
    N = ts.num_individuals
    optima = np.array(ts.metadata['SLiM']['user_metadata']['indsOptimum'])
    mut_effect_lists = []
    for site in ts.sites():
        effect_site = [0]
        for mut in site.mutations:
            effect_site.append(
                mut.metadata['mutation_list'][0]['selection_coeff'])
        mut_effect_lists.append(effect_site)

    # Observed extent of local adaptation (local-foreign contrast (LF))
    (w_local, w_foreign) = LF_fitness(ind_x, ind_y, ind_z,
                                      optima, dist_mate, sigma_w)
    LF_cline = w_local - w_foreign
    mean_LF = np.mean(LF_cline)

    shuffle_replicates = 1
    LF_shuffle = []
    index_site = 0
    # traverse sites with mutations
    for v in ts.variants():
        # Skip calculation and directly append mean_LF as LF_shuffle
        # By definition neutral mutations should not affect phenotypes at all
        if sum(mut_effect_lists[index_site]) == 0:
            index_gt = 1
            while index_gt < len(mut_effect_lists[index_site]):
                LF_shuffle.append(mean_LF)
                index_gt += 1
        else:
            # Calculate the effect of each mutation at the focal site
            index_gt = 1
            # traverse all derived alleles at the site
            while index_gt < len(mut_effect_lists[index_site]):
                # for multiallelic sites, hide other alleles at the site except the focal allele,
                # leave only one allele at each time
                focal_gt = np.array(v.genotypes)
                focal_gt[focal_gt != index_gt] = 0
                # effect of the focal mutation in each individual
                effect_mut_genomes = np.array(mut_effect_lists[index_site])[
                    focal_gt]
                # Add up the two genomes of each individual
                effect_mut_ind = (effect_mut_genomes[range(0, 2 * N - 1, 2)] +
                                  effect_mut_genomes[range(1, 2 * N, 2)])
                phenotype_without_mut = ind_z - effect_mut_ind
                # Shuffle the distribution of the individuals
                # (without changing the observed heterozygosity)
                r = 0
                LF_shuffle_mut = np.zeros(shuffle_replicates)
                while r < shuffle_replicates:
                    effect_mut_ind_shuffle = list(effect_mut_ind)  # shallow copy
                    np.random.shuffle(effect_mut_ind_shuffle)
                    # Phenotype with shuffled mut =
                    # phenotypes of each ind - effect of the focal mutation in each ind
                    # + effect of the focal mutation in each ind after shuffling
                    phenotype_shuffle = phenotype_without_mut + effect_mut_ind_shuffle
                    (w_local, w_foreign) = LF_fitness(ind_x, ind_y,
                                                      phenotype_shuffle,
                                                      optima, dist_mate,
                                                      sigma_w)
                    LF_shuffle_mut[r] = np.mean(w_local - w_foreign)
                    r += 1
                LF_shuffle.append(np.mean(LF_shuffle_mut))
                index_gt += 1
        index_site += 1
        if index_site % 10000 == 0:
            print(f"{index_site} sites processed")
    delta_LF_mut = mean_LF - LF_shuffle
    return(mean_LF, delta_LF_mut)

############################# program #########################################
# Values
sigma_w = 0.4
dist_mate = 0.12
history = 100000 # number of generations before the focal model

path_file_basename = args.input
# path_file_basename = "/home/tianlin/Documents/github/data/slim_data/glacial_history/historical_optimum0_timeSeries/M2b_highPoly_lowMig_clineMap/deleted/Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.03_mateD0.12_seed4444751029864420980_tick"
# path_file_basename = "/home/tianlin/Documents/github/data/slim_data/glacial_history/historical_optimum0_timeSeries/M2b_highPoly_lowMig_clineMap/Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.03_mateD0.12_seed1605901219043795304_tick"
model_name = "_".join(path_file_basename.split("/")[-1].split("_")[0:-2])
seed = path_file_basename.split("_")[-2]
# figPath = "/home/tianlin/Documents/github/data/tskit_data/figure/20240605/"
figPath = "/home/tianlin/Documents/github/data/tskit_data/figure/20240809/one_run/"
# outBasePath = "/home/tianlin/Documents/github/data/tskit_data/output/table/historical_optimum0_timeSeries/"
outBasePath = "/home/tianlin/Documents/github/data/tskit_data/output/table/gradualEnvChange_timeSeries/"
outPath = outBasePath + model_name + "_" + seed + "/"
if not os.path.exists(outPath):
    os.makedirs(outPath)
if not os.path.exists(figPath):
    os.makedirs(figPath)

mean_lf_list = []
# times = [50, 100, 200, 300, 400, 800, 1000, 2000, 4000, 10000,
#          20000, 40000, 100000]
# times = [50, 100, 200, 300, 400, 800, 1000, 2000, 4000, 10000]
# 18 samples
# times = np.array([50, 100, 200, 300, 400, 800,
#                   1000, 2000, 4000, 6000, 8000,
#                   10000, 20000, 30000, 40000, 60000, 80000, 100000])
#22 samples
# times = np.array([0, 100, 200, 400, 600, 800,
#                   1000, 1200, 1400, 1600, 1800,
#                   2000, 4000, 6000, 8000,
#                   10000, 20000, 30000, 40000, 60000, 80000, 100000])

times = np.array([0, 400, 800, 1000, 4000, 10000,
                  20000, 40000, 60000, 80000, 100000])
# times = np.array([10000])
# times = np.array([1000, 4000, 20000, 40000, 60000, 80000, 100000])
LFs = []
gea_goals = []
fileid = 1
corPE = []
for i in times:
    path_file_name = path_file_basename + str(i + 100000) + ".trees"
    # Delete seed and tick information for making a directory
    # Tree-sequence file from SLiM
    ts = tskit.load(path_file_name)
    current_tick = ts.metadata['SLiM']['tick']
    # x,y: coordinates, z: phenotypes
    ind_x, ind_y, phenotypes = zip(*ts.individuals_location)
    # Genomic position of sites
    pos = ts.sites_position
    # Genomic position of mutations
    pos_by_mut = ts.sites_position[ts.mutations_site]
    # Maximum tick
    max_tick = ts.metadata["SLiM"]["tick"]
    # Number of diploid individuals
    N = ts.num_individuals
    # Environmental optima
    optima = np.array(ts.metadata['SLiM']['user_metadata']['indsOptimum'])
    # mapValues = np.array(ts.metadata['SLiM']['user_metadata']['mapValues'])
    # Correlation between phenotypes and optima
    corPE.append(stats.kendalltau(phenotypes, optima)[0])
    # Temporary use
    (w_local, w_foreign) = LF_fitness(ind_x, ind_y, phenotypes,
                                      optima, dist_mate, sigma_w)
    LF_cline = w_local - w_foreign
    mean_LF = np.mean(LF_cline)
    LFs.append(mean_LF)

    # The following part will take a long time:
    # Skip recapitation and adding neutral mutation steps
    # Phenotypic effect of each mutation
    mut_effect = []
    for mut in ts.mutations():
        mut_effect.append(mut.metadata['mutation_list'][0]['selection_coeff'])
    mut_effect = np.array(mut_effect)
    # Age of each mutation
    age = ts.mutations_time
    # Frequency of each mutation
    freq = []
    num_samples = ts.num_samples
    for v in ts.variants():
        for allele in np.arange(1, v.num_alleles):
            focal_freq = np.count_nonzero(v.genotypes == allele) / num_samples
            freq.append(focal_freq)
    # mean LF and LF_mut of each mutation: takes a long time
    mean_lf, delta_LF_mut = calculate_lfmut(ts, dist_mate, sigma_w)
    mean_lf_list.append(mean_lf)
    # Genomic position of mutations
    pos_by_mut = ts.sites_position[ts.mutations_site]
    # Save a table for each time point
    myTable = np.array([pos_by_mut,
                        age,
                        freq,
                        mut_effect,
                        delta_LF_mut])
    np.savetxt(outPath + model_name + "_" + seed + "_tick" + str(current_tick)
               + "_functionalMut_table.txt",
               myTable,
               header="pos, allele_age, allele_freq, effect_size, lf_mut")
    print(str(fileid) + " files processed")
    # Minimum number of alleles for explaining at least 80% positive LF
    expected_explained_proportion = 0.8
    positive_lf_mut = delta_LF_mut[delta_LF_mut > 0]
    sorted_lf_mut = np.array(list(reversed(sorted(positive_lf_mut))))
    positiveTotal = sum(sorted_lf_mut)
    cumulative_lf = [0] + [sum(sorted_lf_mut[0:k + 1]) / positiveTotal
                           for k in range(len(sorted_lf_mut))]
    # How many alleles do we need to account for 80% of current local adaptation?
    gea_goal = sum(np.array(cumulative_lf) < expected_explained_proportion)
    gea_goals.append(gea_goal)

# cor(phe, env) ~ time
plt.figure(figsize=(12, 5))
plt.plot(times, corPE,
         marker="o", markersize=5, markerfacecolor="white",
         color="black")
plt.xlabel("Generations after the environmental change")
plt.ylabel("Correlation between phenotypes and environmental optima "
           "\n(Kendall's tau)")
plt.savefig(figPath + model_name + seed + "_corPE_through_time_100000.png",
            dpi=300)
plt.close()

# LF ~ time
plt.figure(figsize=(12,5))
plt.plot(times, LFs,
         marker="o", markersize=5, markerfacecolor="white",
         color="black")
plt.xlabel("Generations after the environmental change")
plt.ylabel("Extent of local adaptation (LF)")
plt.savefig(figPath + model_name + seed + "_LF_through_time_100000.png",
            dpi=300)
plt.close()

# # gea_goal ~ time
# plt.figure(figsize=(12,5))
# plt.plot(times, gea_goals,
#          marker="o", markersize=5, markerfacecolor="white",
#          color="black")
# plt.xlabel("Generations after the environmental change")
# plt.ylabel("Minimum number of alleles for explaining 80% local adaptation")
# plt.savefig(figPath + model_name + seed + "_geaGoal_through_time_100000.png",
#             dpi=300)
# plt.close()


