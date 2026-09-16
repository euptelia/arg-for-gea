"""
Test the effects of FPR calculation methods
2026.09.15
tianlin.duan42@gmail.com
"""
############################## modules ####################################
import numpy as np
import glob #for loading files
import pandas as pd #dataframe
import os #mkdir
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches #manually make legends

############################# options #############################
import argparse
parser = argparse.ArgumentParser()
# Path should end by "/"
parser.add_argument('-i', '--input',
                    help='Path to input tables',
                    type=str)
# #parser.add_argument('-n', '--name',
#                     help='Short model name',
#                     type=str)
# #parser.add_argument('-p', '--plot',
#                     help='1 for generating plots and 0 for not plotting',
#                     type=int, default=1)
args = parser.parse_args()

############################# program #########################################
# Values
# sigma_w = 0.4
# dist_mate = 0.15
num_runs = 200
inPath = args.input
# inPath = "/media/anadem/PortableSSD/arg4gea_data/tskit_data/tick110000/Continuous_nonWF_M2a_glacialHistoryOptimum0_clineMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.03_mateD0.12_K6000_r1.0e-07/tick110000/"

simName = inPath.split("/")[-3] #Check this before use
#print(simName)

#Short title: Hard coded. Check this before using!
if "sigmaD0.06_mateD0.15" in simName:
    migName = "HighMig"
elif "sigmaD0.03_mateD0.12" in simName:
    migName = "LowMig"
else:
    migName = ""
if "_clineMap_" in simName:
    mapName = "Cline"
elif "_patchyMap_" in simName:
    mapName = "Patchy"
else:
    mapName = ""
if "_sigmaM0.01_" in simName:
    mutName = "highPoly"
elif "_sigmaM0.1_" in simName:
    mutName = "lowPoly"
else:
    mutName = ""
demoName_ori = simName.split("_")[2]
# name_change = {"M2a":"M1a", "M2b":"M1b", "M3a":"M2a", "M3b":"M2b"}
name_change = {"M0a":"Neutral constant","M0b":"Neutral expansion",
               "M2a":"Single change", "M2b":"Single expansion",
               "M3a":"Recurrent changes", "M3b":"Recurrent expansions"}
demoName = name_change[demoName_ori]
shortName = ",".join([demoName, migName, mapName])
name_change2 = {"M0a":"MNeuCon","M0b":"MNeuExp",
                "M2a":"MSelCon", "M2b":"MSelExp",
                "M3a":"MRecCon", "M3b":"MRecExp"}
demoName2 = name_change2[demoName_ori]
# shortName2 = "_".join([demoName2, mutName, migName, mapName])

# figPath = ("/home/anadem/github/data/tskit_data/figure/multiRuns/" +
#            simName + "/" + str(num_runs) + "runs_" +
#            inPath.split("/")[-2]+"/")
figPath = ("/home/anadem/github/data/tskit_data/figure/second_revision/fpr/")
if not os.path.exists(figPath):
    os.makedirs(figPath)

outPath_fpr = ("/home/anadem/github/data/tskit_data/stats/second_revision/fpr_neutralChrom/singleModel/")
if not os.path.exists(outPath_fpr):
    os.makedirs(outPath_fpr)

# Check this before use!
model_name = "_".join(inPath.split("/")[-3:-1] + [str(num_runs)+"runs"])
tick = int(inPath.split("/")[-2].split("tick")[-1])
event_age = tick - 100000
past_event_ages = [tick-i for i in list(range(0, 100000, 10000))]

# Import data
# Shape of the table:
# dataTable[0]     mutation positions
# dataTable[1]     age
# dataTable[2]     freq
# dataTable[3]     mut_effect
# dataTable[4]     delta_LF_mut
# dataTable[5]     cor_GE[0]: Kandell's tau
# dataTable[6]     cor_GE[1]: p-value
# dataTable[7]     relative to the sum of positive delta_LF_mut (delta_LF_mut/sum of positive LF_mut )
# dataTable[8]     relative to sum of delta_LF_mut (delta_LF_mut/sum of LF_mut )
# dataTable[9]     temperary run ID 0,1,2...(number of files-1)
# dataTable[10]    rank of p-values
df = pd.DataFrame()
# Load results of multiple runs from file
fileList = glob.glob(inPath + "*.txt")
#print(inPath)
#print(fileList)
run = 0
for f in fileList:
    # df_focal = pd.read_csv(f, sep='\t', header=0)
    df_focal = pd.read_csv(f, sep="\t", header=0)
    # Add relative delta_LF_mut,temporary run ID, and rank of p-values
    lf_sum_positive = sum(df_focal.loc[df_focal["delta_LF_mut"]>0,"delta_LF_mut"])
    lf_sum_negative = sum(df_focal.loc[df_focal["delta_LF_mut"]<0, "delta_LF_mut"])
    lf_sum = sum(df_focal["delta_LF_mut"])
    df_focal["relative_positive_lf"] = (df_focal["delta_LF_mut"] /
                                        lf_sum_positive)
    df_focal["relative_negative_lf"] = (df_focal["delta_LF_mut"] /
                                        lf_sum_negative)
    df_focal["relative_lf"] = df_focal["delta_LF_mut"] / lf_sum
    df_focal["run_id"] = np.full(shape=(df_focal.shape[0], 1),
                                         fill_value=run)
    # df_focal["p_rank"] = df_focal["p"].argsort().argsort()
    df_focal["p_rank"] = df_focal["p"].rank(ascending=True, pct=True)
    df = pd.concat([df, df_focal], axis=0)
    run += 1
print(str(run) + " files have been loaded.")
# Columns: "pos", "age", "freq", "mut_effect","delta_LF_mut",
#          "tau", "p", "relative_positive_lf", "relative_negative_lf", "relative_lf",
#          "run_id", "p_rank"
del df_focal



#### Test 1: FPR on the neutral chromosome: pos > 7500000
#FPR ~ Allele age (equal time intervals) in neutral alleles
maf_filter = 0.05
# df_plot = df[(df["mut_effect"] == 0) &
#                 (df["freq"] != 0) &
#                 (df["freq"] != 1)]
#Separate neutral alleles with MAF > maf_filter on the neutral chromosome
df_plot = df[(df["mut_effect"] == 0) &
                (df["freq"] >= maf_filter) &
                (df["freq"] <= (1-maf_filter)) &
             (df["pos"] > 7500000)]



# False negative rate for NEUTRAL mutations among AGE categories: Equal intervals
num_cat_age = 22
# num_cat_age = 55
# p_threshold = 0.0001
p_threshold = 0.0000000001
max_age = tick
cat_width_age = max_age/num_cat_age
age_boundaries = np.append(np.arange(0, max_age, cat_width_age), max_age)
FPR_neutral_byAge_equalWidth = []
sample_size_age_equalWidth = []
p_median_byAge_equalWidth = []
p_sd_byAge_equalWidth = []
tau_absMedian_byAge_equalWidth = []
tau_absSd_byAge_equalWidth = []
for i in range(num_cat_age):
    p_category = df_plot["p"][(df_plot["age"] > age_boundaries[i]) &
                                 (df_plot["age"] <= age_boundaries[i+1])]
    tau_category = df_plot["tau"][(df_plot["age"] > age_boundaries[i]) &
                                 (df_plot["age"] <= age_boundaries[i+1])]
    sample_size_age_equalWidth.append(len(p_category))
    p_median_byAge_equalWidth.append(np.nanmedian(p_category))
    p_sd_byAge_equalWidth.append(np.nanstd(p_category))
    tau_absMedian_byAge_equalWidth.append(np.nanmedian(abs(tau_category)))
    tau_absSd_byAge_equalWidth.append(np.nanstd(abs(tau_category)))
    # Neutral mutations have no phenotypic effect and therefore no positive
    mut_in_cat = len(p_category)
    obsP_neutral = p_category < p_threshold
    # All positives are false positives
    FP = sum(obsP_neutral)
    FPR_neutral_byAge_equalWidth.append(FP/mut_in_cat)


# Equal intervals
age_fig_size = (7,5) # 22 bins
# age_fig_size = (8,5) # 40 bins
# age_fig_size = (15,5) # 100 bins
# age_fig_size = (9,5) # 55 bins
label_font = 16
tick_font = 16
x = age_boundaries[0:-1] + cat_width_age/2
y = FPR_neutral_byAge_equalWidth
# Neutral alleles: FPR ～ Allele age, equal intervals
plt.figure(figsize=age_fig_size)
# plt.figure(figsize=(15,5)) # 100 bins
plt.plot(x,
         y,
         color="grey",
         marker = "o")
plt.xlabel("Allele age (thousand ticks)", fontsize=label_font)
plt.ylabel("False positive rate of neutral alleles \n (FP/(FP+TN))",
           fontsize=label_font)
plt.xticks(ticks=age_boundaries,
           labels=[str(int(i/1000)) for i in age_boundaries],
           rotation=90)
plt.tick_params(axis='both', which='major', labelsize=tick_font)
plt.title("_".join([demoName, mutName, migName, mapName]), fontsize=label_font)
plt.tight_layout()
plt.savefig(figPath+model_name +
            str(num_cat_age)+"bins"+"_FPR_vs_age_neutralAllele_GEAp" + str(p_threshold) +
            "_maf" + str(maf_filter) +
            "_equalInterval.png",
            dpi=350)
plt.close()

#Save the FPR~age table for combined figures
out_path_file = (outPath_fpr + model_name +
                 "_p" + str(p_threshold) +
                 "_maf"+str(maf_filter) +
                 "_cat" + str(num_cat_age)+ "_fpr.tab")
header = "\t".join(["age", "fpr", "tick_text"]) + "\n"
with open(out_path_file, "w") as fout:
    fout.write(header)
    for i in range(len(x)):
        outLine = "\t".join([str(x[i]), str(y[i]), str(age_boundaries[i])]) + "\t" +"\n"
        fout.write(outLine)
fout.close()
print("FPR table has been saved.")