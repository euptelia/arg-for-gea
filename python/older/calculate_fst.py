"""
For all .trees files in a folder,
calculate:
    population fst (Weir and Cockman)
    LF (average contrast between local and foreign mean fitness)
Save results in a table
Typically needs about 16 Gb memory
tianlin.duan42@gmail.com
2024.08.24
Last modified
2025.03.18
"""
############################# modules #########################################
import msprime
import tskit
# import pyslim # for recapitation
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.spatial import distance_matrix
import numpy as np
# import pandas
import random
import time
import allel # for allel.weir_cockerham_fst()
import os # mkdir
# import tracemalloc check memory usage
import glob #for loading files

############################# options #############################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inPath',
                    help='Absolute path to the folder with input .trees files',
                    type=str)
parser.add_argument('-o', '--outPath',
                    help='Absolute path for an output table',
                    type=str)
args = parser.parse_args()

############################# program #########################################
#Values
sigma_w = 0.4
dist_mate = 0.15

inPath = args.inPath
outPath = args.outPath
fileList = glob.glob(inPath + "*.trees")
model_name = fileList[0].split(".trees")[0]
if not os.path.exists(outPath):
    os.makedirs(outPath)
output = outPath + "/" + model_name +"fst.tab"
fout = open(output, "w")
for file in fileList:
    # print("Fst calculation started")
    # print(time.ctime())
    file_name = file.split("/")[-1]
    t1 = time.time()
    # Tree-sequence file from SLiM
    ts = tskit.load(file)
    # Genomic position of sites
    pos = ts.sites_position
    # Genomic position of mutations
    pos_by_mut = ts.sites_position[ts.mutations_site]
    # Maximum tick
    max_tick = ts.metadata["SLiM"]["tick"]
    # Number of diploid individuals
    N = ts.num_individuals
    #Location (x,y) and phenotype (z)
    ind_x, ind_y, ind_z = zip(*ts.individuals_location)
    # Environmental optima
    optima = np.array(ts.metadata['SLiM']['user_metadata']['indsOptimum'])
    mapValues = np.array(ts.metadata['SLiM']['user_metadata']['mapValues'])

    # Add neutral mutations
    mutation_seed = random.randint(1, 2**31)
    mut_model1 = msprime.SLiMMutationModel(type=1)
    mts = msprime.sim_mutations(ts,
                                rate=5e-8,
                                random_seed=mutation_seed,
                                model=mut_model1,
                                keep=True) #keep the existing mutations
    del ts

    # Mutation effect in a list of lists,
    # corresponding to phenotypic effect of allele 0/1/2... at each site
    mut_effect_lists = []
    for site in mts.sites():
        effect_site = [0]
        for mut in site.mutations:
            effect_site.append(mut.metadata['mutation_list'][0]['selection_coeff'])
        mut_effect_lists.append(effect_site)

    # A collection of individual IDs for each subpopulation
    inds_subpop = {}
    num_row = 10
    num_col = 10
    map_width = 1.0
    num_sample_pop = num_row * num_col
    sample_pop_ids = np.full(fill_value=np.nan,
                             shape=mts.num_individuals)
    i = 0
    for ind in mts.individuals():
        x = int(ind.location[0]/(map_width/num_row))
        y = int(ind.location[1]/(map_width/num_col))
        # Push the ones on the upper boundaries back
        x = x if x < num_row else (num_row-1)
        y = y if y < num_col else (num_col-1)
        sample_pop_ids[i] = x+y*10
        i += 1

    # A vector of allele frequencies of all mutations for each sample population
    freq_sample_pop = np.full(fill_value=np.nan,
                              shape=(mts.num_mutations, num_sample_pop))
    for focal_pop_id in np.arange(num_sample_pop):
        focal_ind_list = np.array(np.where(sample_pop_ids == focal_pop_id)[0])
        focal_genome_list = np.stack((focal_ind_list*2,
                                      focal_ind_list*2+1)).ravel('F')
        inds_subpop[focal_pop_id] = focal_ind_list

    # Fst
    g = mts.genotype_matrix().reshape(mts.num_sites, mts.num_individuals, 2)
    subpops = [inds_subpop[i] for i in range(100)]
    a, b, c = allel.weir_cockerham_fst(g, subpops)
    # fst_pop = a / (a + b + c)
    fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))
    fout.write(file_name + "\t" + str(round(fst, 4)) + "\n")
    print("Fst:" + str(round(fst, 4)))
    del g, a, b, c, subpops
    print("Fst calculation finished")
    print(time.ctime())
    timer = time.time() - t1
    print(str(round((timer/60),2) )+ " minutes used for calculation.")
fout.close()
