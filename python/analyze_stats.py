"""
Collect Fst and LF by runs and calculate mean and se
#2025.06.09
"""
############################# modules #########################################
import numpy as np
import pandas as pd
import glob
from scipy import stats
############################# functions #######################################
#Short title
def shortenNames(oldName):
    #Replace path/file names with a short model name
    if "sigmaD0.06_mateD0.15" in oldName:
        migName = "HighMig"
    elif "sigmaD0.03_mateD0.12" in oldName:
        migName = "LowMig"
    else:
        migName = ""
    if "_sigmaM0.01_" in oldName:
        polyName = "HighPoly"
    elif "_sigmaM0.1_" in oldName:
        polyName = "LowPoly"
    else:
        polyName = ""
    if "_clineMap_" in oldName:
        mapName = "Cline"
    elif "_patchyMap_" in oldName:
        mapName = "Patchy"
    else:
        mapName = ""
    demoName = oldName.split("/")[-1].split("_")[2]
    shortName = "_".join([demoName, polyName, migName, mapName])
    return shortName
############################# program #########################################
#Fst
#High Migration
fileList_highMig = glob.glob("/home/anadem/github/data/tskit_data/stats/fst/Continuous_nonWF_*sigmaD0.06_mateD0.15*_Fst.tab")
# fileList_highMig = glob.glob("/home/anadem/github/data/tskit_data/wrong_fitness_function/stats/fst/Continuous_nonWF_*sigmaD0.06_mateD0.15*_Fst.tab")
df_list_highMig = [pd.read_csv(fileName, header=None) for fileName in fileList_highMig]
df_fst_highMig = pd.concat(df_list_highMig)
np.mean(df_fst_highMig)
np.std(df_fst_highMig)
stats.sem(df_fst_highMig)

#low Migration
fileList_lowMig = glob.glob("/home/anadem/github/data/tskit_data/stats/fst/Continuous_nonWF_*sigmaD0.03_mateD0.12*_Fst.tab")
# fileList_lowMig = glob.glob("/home/anadem/github/data/tskit_data/wrong_fitness_function/stats/fst/Continuous_nonWF_*sigmaD0.03_mateD0.12*_Fst.tab")
df_list_lowMig = [pd.read_csv(fileName, header=None) for fileName in fileList_lowMig]
df_fst_lowMig = pd.concat(df_list_lowMig)
np.mean(df_fst_lowMig)
np.std(df_fst_lowMig)
stats.sem(df_fst_lowMig)

#LF
fileList = glob.glob("/home/anadem/github/data/tskit_data/stats/lf/Continuous_nonWF_*LF.tab")
df_list_lf = [pd.read_csv(fileName, header=None) for fileName in fileList]

for df, name in zip(df_list_lf, fileList):
    df["modelName"] = shortenNames(name)
    print(stats.ttest_1samp(df[0], 0,alternative="greater"))
df_lf = pd.concat(df_list_lf)

lf_mean=df_lf.groupby(["modelName"])[0].mean()
lf_std=df_lf.groupby(["modelName"])[0].std()
lf_summary = pd.concat([lf_mean,lf_std], axis=1)
lf_summary.columns=["lf_mean", "lf_sd"]

# #Save a summary table
# lf_summary.round(4).to_csv("/home/anadem/github/data/tskit_data/stats/summary/lf_summary.tab", sep="\t")


df_lf_highMig = df_lf.loc[df_lf["modelName"].str.contains("_HighMig_"), 0]
df_lf_lowMig = df_lf.loc[df_lf["modelName"].str.contains("_LowMig_"), 0]
df_lf_patchy = df_lf.loc[df_lf["modelName"].str.contains("_Patchy"), 0]
df_lf_cline = df_lf.loc[df_lf["modelName"].str.contains("_Cline"), 0]
df_lf_highPoly = df_lf.loc[df_lf["modelName"].str.contains("_HighPoly"), 0]
df_lf_lowPoly = df_lf.loc[df_lf["modelName"].str.contains("_LowPoly"), 0]

stats.levene(df_lf_highMig, df_lf_lowMig)
stats.ttest_ind(df_lf_highMig, df_lf_lowMig, equal_var=True)
print(df_lf_highMig.mean())
print(df_lf_lowMig.mean())
print(df_lf_highMig.std())
print(df_lf_lowMig.std())

stats.levene(df_lf_highPoly, df_lf_lowPoly)
stats.ttest_ind(df_lf_highPoly, df_lf_lowPoly, equal_var=True)
print(df_lf_highPoly.mean())
print(df_lf_lowPoly.mean())
print(df_lf_highPoly.std())
print(df_lf_lowPoly.std())

stats.levene(df_lf_cline, df_lf_patchy)
stats.ttest_ind(df_lf_cline, df_lf_patchy, equal_var=False)
print(df_lf_patchy.mean())
print(df_lf_cline.mean())
print(df_lf_patchy.std())
print(df_lf_cline.std())




