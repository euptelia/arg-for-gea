############################# modules #########################################
import msprime
import pandas
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
from collections import Counter
import pandas #dataframe

############################# program #########################################
# Values
# model_name = "Continuous_nonWF_M2b_neutralHistory_mu1.0e-09_sigmaM0.1_sigmaW0.4_tick20000_100runs"
# model_name = "Continuous_nonWF_M2b_glacialHistory_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.1_mateD0.2_tick1000_100runs"
model_name = "Continuous_nonWF_M2b_glacialHistory_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.1_mateD0.2_tick10000_100runs"
# sigma_w = 0.4
sigma_w = 1.0
# sigma_w = 2.0
dist_mate = 0.1

# inPath = "/home/tianlin/Documents/github/data/tskit_data/output/table/m2b_mu1e-7_m0.1_w0.4/"
inPath = "/home/tianlin/Documents/github/data/tskit_data/output/table/Continuous_nonWF_M2b_glacialHistory_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.1_mateD0.2/tick110000/"
figPath = "/home/tianlin/Documents/github/data/tskit_data/figure/20240418/100runs_tick10000/"
outPath = "/home/tianlin/Documents/github/data/tskit_data/output/"

# Shape of the table:
# dataTable[0]     mutation ID
# dataTable[1]     age,
# dataTable[2]     freq,
# dataTable[3]     mut_effect,
# dataTable[4]     delta_LF_mut,
# dataTable[5]     cor_GE[0]: Kandell's tau
# dataTable[6]     cor_GE[1]: p-value
# dataTable[7]     relative to the sum of positive delta_LF_mut (delta_LF_mut/sum of positive LF_mut )
# dataTable[8]     relative to sum of delta_LF_mut (delta_LF_mut/sum of LF_mut )
# dataTable[9]     temperary run ID 0,1,2...(number of files-1)
# dataTable[10]    rank of p-values
dataTable = np.empty((11, 0))
# Load results of multiple runs from file
fileList = glob.glob(inPath + "*.txt")
run = 0
for f in fileList:
    focalTable = np.loadtxt(f)
    # print(np.shape(focalTable))
    # Add relative positive delta_LF_mut (delta_LF_mut/sum of positive LF_mut)
    relative_positive_lf_mut = focalTable[4]/sum(focalTable[4][focalTable[4]>0])
    relative_lf_mut = focalTable[4] / sum(focalTable[4])
    p_rank = focalTable[6].argsort().argsort()
    # print(len(focalTable[0]))
    # Add relative delta_LF_mut ,temporary run ID, and rank of p-values
    focalTable = np.concatenate((focalTable,
                                 [relative_positive_lf_mut],
                                 [relative_lf_mut],
                                 np.full(shape=(1, len(focalTable[0])),
                                         fill_value=run),
                                 [p_rank]),
                                axis=0)
    dataTable = np.concatenate((dataTable, focalTable), axis=1)
    run += 1



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

# Change to pandas dadaframe
dataTable_t = dataTable.transpose()
df = pandas.DataFrame(data=dataTable_t,
                      columns=["id", "age", "freq",
                               "mut_effect","delta_LF_mut","tau",
                               "p", "relative_positive_lf", "relative_lf",
                               "run_id", "p_rank"])

# Plots
# Relative LF_mut ~ cor_GE p-values RANKS, colored by age rank bins
p_rank = df["p_rank"]
age_rank_bin = (dataTable[1].argsort().argsort()/10000).astype(int)
plt.scatter(p_rank, dataTable[7],
         marker="o",
         c=age_rank_bin, alpha=0.05)
plt.xlabel("GEA p-values ranks")
plt.ylabel("Relative $LF_{mut}$")
plt.colorbar().set_label("Age rank bins \n (10000 alleles per bin)")
plt.tight_layout()
plt.savefig(figPath+model_name+"_corGEpvalueRank_vs_relativeLFmut_ageBinColor2.png",
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


#### Cumulative plot of LF_mut ####
expected_explained_proportion = 0.8
num_mut_runs = Counter(dataTable[9]).values()

cumulative_lf_list = []
gea_goal = []
explained = []
for i in range(100):
    # print(i)
    focal_relative_lf_mut = dataTable[7][dataTable[9] == i]
    relative_positive_lf_mut = focal_relative_lf_mut[focal_relative_lf_mut > 0]
    sorted_lfmut = np.array(list(reversed(sorted(relative_positive_lf_mut))))
    # Trim the tails with near-0 contributions to ensure all runs have the same
    # number of mutations
    #sorted_lfmut = sorted_lfmut[0:minimum_muts-1]
    # positiveTotal = sum(sorted_lfmut)
    focal_cumulative_lf = [0] + [sum(sorted_lfmut[0:k+1])
                                 for k in range(len(sorted_lfmut))]
    # How many alleles do we need to account for 80% of current local adaptation?
    focal_gea_goal = sum(np.array(focal_cumulative_lf) < expected_explained_proportion)
    focal_explained = focal_cumulative_lf[focal_gea_goal]
    cumulative_lf_list.append(focal_cumulative_lf)
    gea_goal.append(focal_gea_goal)
    explained.append(focal_explained)

# Convert the "list of lists" to a numpy array by adding np.nan values to the end
max_num_positive_mut = len(max(cumulative_lf_list, key=len))
cumulative_lf = np.full(shape=(run, max_num_positive_mut), fill_value=np.nan)
for i in range(run):
    n = len(cumulative_lf_list[i])
    cumulative_lf[i][0:n] = cumulative_lf_list[i]

# mean_cumulative_lf = np.mean(cumulative_lf, axis=0)
# median_cumulative_lf = np.median(cumulative_lf, axis=0)
# mean_trim_len = np.count_nonzero(~np.isnan(mean_cumulative_lf))
# median_trim_len = np.count_nonzero(~np.isnan(median_cumulative_lf))

# Calculate mean and meadian and ignore np.nan values
mean_cumulative_lf = np.nanmean(cumulative_lf, axis=0)
median_cumulative_lf = np.nanmedian(cumulative_lf, axis=0)
mean_gea_goal = int(np.mean(gea_goal))
mean_explained = mean_cumulative_lf[mean_gea_goal]
median_gea_goal = int(np.mean(gea_goal))
median_explained = median_cumulative_lf[median_gea_goal]

# Plot cumulative_lf: Mean with range
for i in range(run):
    y = np.array(cumulative_lf_list[i])
    x = range(len(y))
    plt.plot(x, y,
             color="lightseagreen", alpha = 0.05)
plt.plot(range(len(mean_cumulative_lf)), mean_cumulative_lf,
             color="black")
plt.xlabel("Mutations sorted in descending order of $LF_{mut}$")
plt.ylabel("Cumulative proportion of positive LF")
    # plt.axvline(x=mean_gea_goal, color="firebrick", linestyle="dotted")
    # plt.axhline(y=mean_explained, color="firebrick", linestyle="dotted")

plt.plot([mean_gea_goal, mean_gea_goal], [0, mean_explained],
         color="firebrick", linestyle="dotted")
plt.plot([0, mean_gea_goal], [mean_explained, mean_explained],
         color="firebrick", linestyle="dotted")
plt.annotate(str(mean_gea_goal), xy=(mean_gea_goal+1, mean_explained),
             xytext=(mean_gea_goal + 1 + len(mean_cumulative_lf)/40, 0),
             color="firebrick")
plt.savefig(figPath + model_name +
            "_meanLFcumulative_positive_100runs_meanAndRange.png",
            dpi=300)
plt.close()

# Plot cumulative_lf: Median with range
for i in range(run):
    y = np.array(cumulative_lf_list[i])
    x = range(len(y))
    plt.plot(x, y,
             color="lightseagreen", alpha = 0.05)
plt.plot(range(len(median_cumulative_lf)), median_cumulative_lf,
             color="black")
plt.xlabel("Mutations sorted in descending order of $LF_{mut}$")
plt.ylabel("Cumulative proportion of positive LF")
plt.plot([median_gea_goal, median_gea_goal], [0, median_explained],
         color="firebrick", linestyle="dotted")
plt.plot([0, median_gea_goal], [median_explained, median_explained],
         color="firebrick", linestyle="dotted")
plt.annotate(str(median_gea_goal), xy=(median_gea_goal+1, median_explained),
             xytext=(median_gea_goal + 1 + len(median_cumulative_lf)/40, 0),
             color="firebrick")
plt.savefig(figPath + model_name +
            "_medianLFcumulative_positive_100runs_medianAndRange.png",
            dpi=300)
plt.close()