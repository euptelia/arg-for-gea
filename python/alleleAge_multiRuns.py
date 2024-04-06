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
import glob #for loading files

############################# program #########################################
# Values
model_name = "Continuous_nonWF_M2b_neutralHistory_mu1.0e-09_sigmaM0.1_sigmaW0.4_tick20000_100runs"
# sigma_w = 0.4
sigma_w = 1.0
# sigma_w = 2.0
dist_mate = 0.1

inPath = "/home/tianlin/Documents/github/data/tskit_data/output/table/m2b_mu1e-7_m0.1_w0.4/"
figPath = "/home/tianlin/Documents/github/data/tskit_data/figure/"
outPath = "/home/tianlin/Documents/github/data/tskit_data/output/"

# Shape of tables:
# dataTable[0]     mutation ID
# dataTable[1]     age,
# dataTable[2]     freq,
# dataTable[3]     mut_effect,
# dataTable[4]     delta_LF_mut,
# dataTable[5]     cor_GE[0]: Kandell's tau
# dataTable[6]     cor_GE[1]: p-value
# dataTable[7]     relative delta_LF_mut (delta_LF_mut/sum of positive LF_mut)
#                  Added by this script
dataTable = np.empty((8, 0))
# Load results of multiple runs from file
fileList = glob.glob(inPath + "*.txt")
for f in fileList:
    focalTable = np.loadtxt(f)
    # print(np.shape(focalTable))
    # Add relative delta_LF_mut (delta_LF_mut/sum of positive LF_mut)
    relative_lf_mut = focalTable[4]/sum(focalTable[4][focalTable[4]>0])
    # print(len(focalTable[0]))
    focalTable = np.concatenate((focalTable,
                                np.full(shape=(1, len(focalTable[0])),
                                        fill_value=relative_lf_mut)),
                                axis=0)
    dataTable = np.concatenate((dataTable, focalTable), axis=1)

# True positive rate, False negative rate, False discovery rate
# True adaptive allele: account for more than LF_min percent of positve LF_mut
# GEA significant: BH-adjusted p-value < 0.05
LF_min = 0.01
num_cat = 10
cat_width = int(100/num_cat) if 100 % num_cat == 0 else 100/num_cat
age_percentile = np.percentile(dataTable[1],
                               np.append(np.arange(0, 100, cat_width),
                                         100))
TPR = []
FPR = []
FDR = []
for i in range(num_cat):
    focal_age = np.logical_and(dataTable[1] > age_percentile[i],
                               dataTable[1] <= age_percentile[i+1])
    dataTable_category = dataTable[:, focal_age]
    # num_allele = len(p_BH_category)
    # proportion_sig.append(sum(p_BH_category < 0.05)/num_allele)
    expP = dataTable_category[7] > LF_min
    obsP = dataTable_category[6] < 0.05
    TP = sum(expP & obsP)
    FP = sum(~expP & obsP)
    FN = sum(expP & ~obsP)
    TN = sum(~expP & ~obsP)
    TPR.append(TP/(TP+FN))
    FPR.append(FP/(TN+FP))
    FDR.append(FP/(TP+FP))



# Plots
# Relative LF_mut ~ cor_GE p-values RANKS, colored by age rank bins
p_rank = dataTable[6].argsort().argsort()
age_rank_bin = (dataTable[1].argsort().argsort()/10000).astype(int)
plt.scatter(p_rank, dataTable[7],
         marker="o",
         c=age_rank_bin, alpha=0.2)
plt.xlabel("GEA p-values ranks")
plt.ylabel("Relative $LF_{mut}$")
plt.colorbar().set_label("Age rank bins \n (10000 alleles per bin)")
plt.tight_layout()
plt.savefig(figPath+model_name+"_corGEpvalueRank_vs_relativeLFmut_ageBinColor.png",
            dpi=300)
plt.close()

# By allele age bins
# Relative LF_mut ~ cor_GE p-values RANKS
num_cat = 20
age_percentile = np.percentile(dataTable[1],
                               np.append(np.arange(0, 100, 100/num_cat),
                                         100))
lf_min = min(dataTable[7])
lf_max = max(dataTable[7])
for i in range(num_cat):
    p_category = dataTable[6][np.logical_and(dataTable[1] > age_percentile[i],
                                          dataTable[1] <= age_percentile[i+1])]
    p_category_rank = p_category.argsort().argsort()
    relative_LF_mut = dataTable[7][
        np.logical_and(dataTable[1] > age_percentile[i],
                       dataTable[1] <= age_percentile[i + 1])
    ]
    num_allele = len(p_category)
    # proportion_sig.append(sum(p_category < 0.05) / num_allele)
    plt.scatter(p_category_rank, relative_LF_mut,
                marker="o",
                c="grey", alpha=0.2)
    plt.ylim(lf_min*1.10, lf_max*1.05)
    plt.xlabel("GEA p-values ranks")
    plt.ylabel("Relative $LF_{mut}$")
    plt.title("Rank bin" + str(i))
    plt.tight_layout()
    plt.savefig(figPath + model_name + "AgeBin" + str(i) +
                "_corGEpvalueRank_vs_relativeLFmut.png",
                dpi=300)
    plt.close()


num_cat = 10
# TPR ～ Allele age
plt.plot(np.arange(100/num_cat, 101, 100/num_cat),
         TPR,
         color="grey")
plt.xlabel("Allele age percentile bins")
# plt.ylabel("Proportion of alleles with BH-adjusted p-value < 0.05")
plt.ylabel("True positive rate (TP/(TP+FN))")
plt.xticks(ticks=np.arange(100/num_cat, 101, 100/num_cat),
           labels=[(str(i*cat_width)+"-"+str((i+1)*cat_width))
                   for i in range(num_cat)])
# plt.savefig(figPath+model_name+"_proportionGEABHp_vs_age.png",
#             dpi=300)
plt.tight_layout()
plt.savefig(figPath+model_name+str(num_cat)+"_minLF"+str(LF_min)+"_"+
            str(num_cat)+"bins"+"_TPR_lfmut1percent_GEABHp0.05_vs_age.png",
            dpi=300)
plt.close()

# FPR ～ Allele age
plt.plot(np.arange(100/num_cat, 101, 100/num_cat),
         FPR,
         color="grey")
plt.xlabel("Allele age percentile bins")
# plt.ylabel("Proportion of alleles with BH-adjusted p-value < 0.05")
plt.ylabel("False positive rate (FP/(FP+TN)))")
plt.xticks(ticks=np.arange(100/num_cat, 101, 100/num_cat),
           labels=[(str(i*cat_width)+"-"+str((i+1)*cat_width))
                   for i in range(num_cat)])
# plt.savefig(figPath+model_name+"_proportionGEABHp_vs_age.png",
#             dpi=300)
plt.tight_layout()
plt.savefig(figPath+model_name+str(num_cat)+"_minLF"+str(LF_min)+"_"+
            str(num_cat)+"bins"+"_FPR_lfmut1percent_GEABHp0.05_vs_age.png",
            dpi=300)
plt.close()


# FDR ～ Allele age
plt.plot(np.arange(100/num_cat, 101, 100/num_cat),
         FDR,
         color="grey")
plt.xlabel("Allele age percentile bins")
# plt.ylabel("Proportion of alleles with BH-adjusted p-value < 0.05")
plt.ylabel("False discovery rate (FP/(FP+TP)))")
plt.xticks(ticks=np.arange(100/num_cat, 101, 100/num_cat),
           labels=[(str(i*cat_width)+"-"+str((i+1)*cat_width))
                   for i in range(num_cat)])
# plt.savefig(figPath+model_name+"_proportionGEABHp_vs_age.png",
#             dpi=300)
plt.tight_layout()
plt.savefig(figPath+model_name+str(num_cat)+"_minLF"+str(LF_min)+"_"+
            str(num_cat)+"bins"+"_FDR_lfmut1percent_GEABHp0.05_vs_age.png",
            dpi=300)
plt.close()