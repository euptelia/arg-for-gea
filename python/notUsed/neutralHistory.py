#!/usr/bin/env python3
"""Modified from SLiM Manual-3, section 17.8
Make neutral coalescent history and add neutral mutations to trees
Generate .tree input for SLiM as a neutral burn-in
！NOTE: check the genome structure of SLiM model before use!
tianlin.duan42@gmail.com
2024.04.05
"""
############################# modules #############################
import msprime
import pyslim
import numpy as np
from time import time
############################# options #############################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-N', '--populationSize',
                    help='population size (number of diploid individuals)',
                    type=int, default=5000)
parser.add_argument('-n', '--samples',
                    help='number of samples',
                    type=int, default=5000)
parser.add_argument('-L', '--seqLen',
                    help='Length of the sequence (genome)',
                    type=int, default=1e7)
parser.add_argument('-r', '--recombinationRate',
                    help='recombination rate',
                    type=float, default=1e-7)
parser.add_argument('-mu', '--mutationRate',
                    help='mutation rate',
                    type=float, default=1e-8)
parser.add_argument('-rep', '--replicates',
                    help='Number of replicate',
                    type=int, default=1)
args = parser.parse_args()

############################# program #############################
N = int(args.populationSize)
n = int(args.samples)
L = int(args.seqLen)
r = float(args.recombinationRate)
mu = float(args.mutationRate)
num_replicates = int(args.replicates)
rng = np.random.RandomState(42)
mut_model = msprime.SLiMMutationModel(type=2)
saveDir = "/home/tianlin/Documents/github/data/slim_data/msprime_history/"
# saveDir = "D:/UBC_office/slim/input/"
t = time()

# Match the recombination landscape in SLiM: Four chromosomes
# [   r   ] 0.5 [   r   ] 0.5 [   r   ] 0.5 [   r   ]
len_chrom = L/4
breaks = [0, len_chrom - 1,
          len_chrom, len_chrom * 2 - 1,
          len_chrom * 2, len_chrom * 3 - 1,
          len_chrom * 3, len_chrom * 4]
rates = [r, 0.5] *  3 + [r]
recomb_map = msprime.RateMap(position=breaks, rate=rates)
print(recomb_map)

# Neutral demographic history
neutral_demography = msprime.Demography()
neutral_demography.add_population(name="p1", initial_size=N)

for i in range(num_replicates):
    (ancestry_seed, mutation_seed) = rng.randint(1, 2**31, size=2)
    # Add neutral coalescent history
    ts = msprime.sim_ancestry(samples=n,
                              sequence_length=L,
                              recombination_rate=recomb_map,
                              random_seed=ancestry_seed,
                              demography=neutral_demography)
    slim_ts = pyslim.annotate(ts, model_type="nonWF", tick=1, stage="early")

    # Match the mutation map in SLiM:
    breaks = []
    len_chrom = int(L / 4)
    len_gene = 5000
    len_intergenic_list = [10000, 10000, 10000]
    # gene are embedded in intergenic regions with different lengths
    for j in range(3):
        len_intergenic = len_intergenic_list[j]
        len_block = len_gene + len_intergenic
        for k in range(len_chrom * j,
                       int(len_chrom * (j + 1) - len_block),
                       len_block):
            breaks.extend([k, k + len_gene])
    mut_map = msprime.RateMap(position=breaks+[L],
                              rate=[mu, 0]*int(len(breaks)/2))
    print(mut_map)
    print(len(mut_map))

    # Add m2 mutations
    mts = msprime.sim_mutations(slim_ts,
                                rate=mut_map,
                                model=mut_model,
                                random_seed=mutation_seed,
                                keep=True)
    mts.dump(saveDir + "neutral_" +
                 "N" + str(N) + "_L" + str(L) + "_r" + str(r) + "_mu" + str(mu) +
                 "_ancSeed" + str(ancestry_seed) + "_mutSeed" + str(mutation_seed) +
                 "_rep" + str(i+1) + "_intergene10kbx3chrom.trees")
timer = (time() - t)/60
print(timer)  # 0.20
