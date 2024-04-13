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

inPath = "/home/tianlin/Documents/github/data/slim_data/"
# inPath = "/home/tianlin/Documents/github/data/slim_data/m2b_mu1e-7_m0.1_w0.4/test/"
figPath = "/home/tianlin/Documents/github/data/tskit_data/figure/"
outPath = "/home/tianlin/Documents/github/data/tskit_data/output/"
# model_name = "Continuous_nonWF_M2a_msprimeHistory_randomMap_mu1.0e-08_r1.0e-07_sigmaM0.01_sigmaW0.4_sigmaD0.02_mateD0.1_seed4448450953634440252_tick100"
model_name = "Continuous_nonWF_M2a_msprimeHistory_randomMap_mu1.0e-08_r1.0e-07_sigmaM0.01_sigmaW0.4_sigmaD0.02_mateD0.1_seed157495413993534526_tick100"
ts = tskit.load(inPath + model_name + ".trees")
