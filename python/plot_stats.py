"""
Analyze/Plot number of individuals, individual age, k80, fst from simulations
tianlin.duan42@gmail.com
2025.03.19
Last modified on
2025.11.10
"""
############################# modules #########################################
import numpy as np
import random
from time import time
import sys # for sys.exit()
import glob #for loading files
from collections import Counter
import pandas as pd#dataframe
import os #mkdir
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import f_oneway
import scipy.stats as stats
import matplotlib.patches as mpatches # manually make legends

############################# program #########################################
#### N individual ####
n_ind = pd.read_table("/home/anadem/github/data/tskit_data/stats/test/n_ind/n_ind_32model.tab",
                    sep="\t")
figPath = "/home/anadem/github/data/tskit_data/figure/multiModels/"
if not os.path.exists(figPath):
    os.makedirs(figPath)
n_ind.mean()
n_ind.std()
n_ind.mean(axis=None)
n_ind.stack().std(axis=None)

# Individual ages
# Individual age of 32 models
inPath = "/home/anadem/github/data/tskit_data/stats/test/ind_age/"
outPath = "/home/anadem/github/data/tskit_data/stats/"
fileList = glob.glob(inPath + "*.tab")
ind_age_dict = {}
for f in fileList:
    df_focal = pd.read_csv(f, sep="\t", header=0)
    focal_name=list(df_focal)[0]
    ind_age_dict[focal_name]=df_focal[focal_name]
ind_age = pd.DataFrame.from_dict(ind_age_dict)
ind_age = ind_age.reindex(sorted(ind_age.columns), axis=1)
del ind_age_dict
# ind_age.mean()
ind_age.stack().mean()
ind_age.median()
# ind_age.quantile(0.75)
# ind_age.quantile(0.8)
ind_age.quantile(0.95)
ind_age.quantile(0.99)
ind_age.quantile(0.999)
ind_age.max()
#Proportion of individuals lives longer than 200 years in percentage
round((ind_age > 200).sum()/(len(ind_age.index)-np.isnan(ind_age).sum())*100, 3)
round((ind_age > 300).sum()/(len(ind_age.index)-np.isnan(ind_age).sum())*100, 3)
round((ind_age > 1000).sum()/(len(ind_age.index)-np.isnan(ind_age).sum())*100, 3)

#Overall
ind_age.stack().quantile(0.999)
(ind_age.stack() > 200).sum()/(len(ind_age.stack())-np.isnan(ind_age.stack()).sum())
(ind_age.stack() > 300).sum()/(len(ind_age.stack())-np.isnan(ind_age.stack()).sum())

#Write the summary table
myTable = pd.DataFrame({"Model": ind_age.columns,
                        "Median": ind_age.median(),
                        "95th_percentile": round(ind_age.quantile(0.95), 3),
                        "Maximum": ind_age.max(),
                        "Percentage_of_individuals_older_than_200_years_old_(%)": round((ind_age > 200).sum() / (len(ind_age.index) - np.isnan(ind_age).sum()) * 100, 3)
                        })
myTable.to_csv(outPath+"TableSX_age_summary.txt", sep='\t', encoding='utf-8', index=False, header=True)

#Log10 boxplot: not accurate as age 0 is excluded
fig, axs = plt.subplots(nrows=1, ncols=1,
                        figsize=(12,6))
# axs.boxplot(ind_age)
np.log10(ind_age).boxplot(notch=True, grid=True, rot=90)
fig.tight_layout()
plt.savefig(figPath+
            "individual_age_boxplot.png",
            dpi=300)
plt.close()

np.log10(ind_age)

# #Violin plot: not used
# # Remove NaN values from the dataset before using matplotlib
# ind_age_noNan = [np.array(np.log10(ind_age[group]))[~np.isnan(np.log10(ind_age[group]))] for group in ind_age]
# plt.figure(figsize=(12,6))
# plt.violinplot(ind_age_noNan)
# plt.xticks(ticks=np.arange(len(ind_age_noNan))+1.125,
#            labels=list(ind_age),
#            rotation=90,
#            ha='center')
# fig.tight_layout()
# plt.savefig(figPath+
#             "individual_age_violinPlot_log10.png",
#             dpi=300)
# plt.close()



#### k80 ####
k80 = pd.read_table("/home/anadem/github/data/tskit_data/stats/k80/k80_32model.tab",
                    sep="\t")
# figPath = "/home/anadem/github/data/tskit_data/figure/multiRuns/"
figPath = "/home/anadem/github/data/tskit_data/figure/multiModels/"
if not os.path.exists(figPath):
    os.makedirs(figPath)

# select columns containing 'lowPoly'
k80_lowPoly = k80.filter(like="LowPoly", axis=1)
k80_highPoly = k80.filter(like="HighPoly", axis=1)
k80_lowMig = k80.filter(like="LowMig", axis=1)
k80_highMig = k80.filter(like="HighMig", axis=1)
k80_lowPoly_lowMig = k80_lowPoly.filter(like="LowMig", axis=1)
k80_lowPoly_highMig = k80_lowPoly.filter(like="HighMig", axis=1)
k80_highPoly_lowMig = k80_highPoly.filter(like="LowMig", axis=1)
k80_highPoly_highMig = k80_highPoly.filter(like="HighMig", axis=1)

#Mean and sd
#Polygenicity
print(k80_lowPoly.mean(axis=None))
print(k80_highPoly.mean(axis=None))
print(k80_lowPoly.stack().std())
print(k80_highPoly.stack().std())
# print(k80_lowPoly.stack().sem())
# print(k80_highPoly.stack().sem())

#expansion in highPoly
# print(k80.filter(like="Con").mean(axis=None))
# print(k80.filter(like="Exp").mean(axis=None))
# print(k80.filter(like="Con").stack().std())
# print(k80.filter(like="Exp").stack().std())
# print(k80.filter(like="Con").stack().sem())
# print(k80.filter(like="Exp").stack().sem())
print(k80_highPoly.filter(like="Con").mean(axis=None))
print(k80_highPoly.filter(like="Exp").mean(axis=None))
print(k80_highPoly.filter(like="Con").stack().std())
print(k80_highPoly.filter(like="Exp").stack().std())
# print(k80_highPoly.filter(like="Con").stack().sem())
# print(k80_highPoly.filter(like="Exp").stack().sem())
#Single expansion
print(k80_highPoly.filter(like="Con").filter(like="Sel").mean(axis=None))
print(k80_highPoly.filter(like="Exp").filter(like="Sel").mean(axis=None))
print(k80_highPoly.filter(like="Con").filter(like="Sel").stack().std())
print(k80_highPoly.filter(like="Exp").filter(like="Sel").stack().std())
print((1-k80_highPoly.filter(like="Exp").filter(like="Sel").mean(axis=None)/k80_highPoly.filter(like="Con").filter(like="Sel").mean(axis=None))*100)
#Recurrent expansion
print(k80_highPoly.filter(like="Con").filter(like="Rec").mean(axis=None))
print(k80_highPoly.filter(like="Exp").filter(like="Rec").mean(axis=None))
print(k80_highPoly.filter(like="Con").filter(like="Rec").stack().std())
print(k80_highPoly.filter(like="Exp").filter(like="Rec").stack().std())
print((1-k80_highPoly.filter(like="Exp").filter(like="Rec").mean(axis=None)/k80_highPoly.filter(like="Con").filter(like="Rec").mean(axis=None))*100)


#Migration
print(k80_lowMig.stack().mean(axis=None))
print(k80_highMig.stack().mean(axis=None))
print(k80_lowMig.stack().std()) #
print(k80_highMig.stack().std()) #
print((1-k80_highMig.stack().mean(axis=None)/k80_lowMig.stack().mean(axis=None))*100)
# print(k80_lowMig.stack().sem()) #
# print(k80_highMig.stack().sem()) #
print(k80_highPoly_lowMig.stack().mean(axis=None))
print(k80_highPoly_highMig.stack().mean(axis=None))
print(k80_highPoly_lowMig.stack().std()) #
print(k80_highPoly_highMig.stack().std()) #
print((1-k80_highPoly_highMig.stack().mean(axis=None)/k80_highPoly_lowMig.stack().mean(axis=None))*100)
# print(k80_lowPoly_lowMig.stack().mean(axis=None))
# print(k80_lowPoly_highMig.stack().mean(axis=None))
# print(k80_lowPoly_lowMig.stack().sem()) #
# print(k80_lowPoly_highMig.stack().sem()) #

#cline/patchy
print(k80.filter(like="Cline").mean(axis=None)) #
print(k80.filter(like="Patchy").mean(axis=None)) #
print(k80.filter(like="Cline").stack().std()) #
print(k80.filter(like="Patchy").stack().std()) #

print(k80_highPoly.filter(like="Cline").mean(axis=None)) #
print(k80_highPoly.filter(like="Patchy").mean(axis=None)) #
print(k80_highPoly.filter(like="Cline").stack().sem()) #
print(k80_highPoly.filter(like="Patchy").stack().sem()) #
print((1-k80_highPoly.filter(like="Cline").mean(axis=None)/k80_highPoly.filter(like="Patchy").mean(axis=None))*100)


#two samplet test
#Demographis history
stats.levene(k80_highPoly.filter(like="Con").stack(),
                k80_highPoly.filter(like="Exp").stack())
stats.ttest_ind(k80_highPoly.filter(like="Con").stack(),
                k80_highPoly.filter(like="Exp").stack(),
                equal_var=False)

stats.ttest_ind(k80_highPoly.filter(like="Con").stack(),
                k80_highPoly.filter(like="Con").filter(like="Sel").stack(),
                equal_var=False)

stats.ttest_ind(k80_highPoly.filter(like="Con").stack(),
                k80_highPoly.filter(like="Con").filter(like="Rec").stack(),
                equal_var=False)

#Migration
stats.levene(k80_highPoly_lowMig.stack(),
                k80_highPoly_highMig.stack())
stats.ttest_ind(k80_highPoly_lowMig.stack(),
                k80_highPoly_highMig.stack(),
                equal_var=False)

#Cline vs. Patchy
stats.levene(k80_highPoly.filter(like="Cline").stack(),
                k80_highPoly.filter(like="Patchy").stack())
stats.ttest_ind(k80_highPoly.filter(like="Cline").stack(),
                k80_highPoly.filter(like="Patchy").stack(),
                equal_var=False)



# #One way ANOVA
# f_oneway(k80.filter(like="Con").stack(),
#          k80.filter(like="Exp").stack())
# f_oneway(k80_lowMig.stack(),
#          k80_highMig.stack())
# f_oneway(k80.filter(like="Cline").stack(),
#          k80.filter(like="Patchy").stack())
# #Only high poly
# f_oneway(k80_highPoly.filter(like="Con").stack(),
#          k80_highPoly.filter(like="Exp").stack())
# f_oneway(k80_highPoly_lowMig.stack(),
#          k80_highPoly_highMig.stack())
# f_oneway(k80_highPoly.filter(like="Cline").stack(),
#          k80_highPoly.filter(like="Patchy").stack())



#High polygenicity
# focal_data1 = k80_highPoly_lowMig
# focal_data2 = k80_highPoly_highMig
# #low migration, high migration
# colors = ["#4a2522","#a68785"]
# patch1 = mpatches.Patch(color=colors[0], label='Low migration')
# patch2 = mpatches.Patch(color=colors[1], label='High migration')
# widths = 0.2
# labels1 = [".".join(i.split("_")[0:2])[1:-3] for i in list(focal_data1)]
# labels2 = [".".join(i.split("_")[0:2])[1:-3] for i in list(focal_data2)]
# labels = [i for j in list(zip(labels1, labels2)) for i in j]

# #Single Plot
# fig, ax = plt.subplots()
# ax.boxplot(focal_data1, notch=True, patch_artist=True, widths=widths,
#            boxprops=dict(facecolor=colors[0], color=colors[0]),
#            capprops=dict(color=colors[0]),
#            whiskerprops=dict(color=colors[0]),
#            flierprops=dict(color=colors[0], markeredgecolor=colors[0]),
#            medianprops=dict(color=colors[0]))
# ax.boxplot(focal_data2, positions=np.arange(len(focal_data2.columns))+1.3,
#            notch=True, patch_artist=True, widths=widths,
#            boxprops=dict(facecolor=colors[1], color=colors[1]),
#            capprops=dict(color=colors[1]),
#            whiskerprops=dict(color=colors[1]),
#            flierprops=dict(color=colors[1], markeredgecolor=colors[1]),
#            medianprops=dict(color=colors[1])
#            )
# focal_min = np.concat([focal_data1, focal_data2], axis=1).min(axis=None)
# focal_max = np.concat([focal_data1, focal_data2], axis=1).max(axis=None)
# focal_range = focal_max-focal_min
# plt.ylim(focal_min - focal_range*0.1,
#         focal_max + focal_range*0.2)
# # ax.get_xticks()
# print(ax.get_xticks())
# plt.xticks(ticks=np.arange(len(focal_data1.columns))+1.125,
#            labels=labels1,
#            # rotation=20,
#            ha='center')
# plt.xlabel(xlabel="Models")
# plt.ylabel(ylabel="Minimum number of alleles for explaining \n"+r"80% allelic contribution to local adaptation ( $\it{K_{80}}$)")
# plt.legend(handles=[patch2, patch1],
#            title="",
#            loc="upper center")
# plt.savefig(figPath+
#             "k80_20models_boxplot_highPoly.png",
#             dpi=300)
# # (left, bottom, right, top),
# plt.tight_layout(rect=(0, 0, 1, 1))
# plt.close()


#Two plots: Low polygenicity and high polygenicity
colors = ["#4a2522","#a68785"]
widths = 0.2
patch1 = mpatches.Patch(color=colors[0], label='Low migration')
patch2 = mpatches.Patch(color=colors[1], label='High migration')

fig, axs = plt.subplots(nrows=2, ncols=1,
                        figsize=(5,9))
#Upper plot: high polygenicity
focal_data1 = k80_highPoly_lowMig
focal_data2 = k80_highPoly_highMig
labels_map = [i.split("_")[3] for i in list(focal_data1)]
labels_model = ["\n\n"+i.split("_")[0] for i in list(focal_data1)][0::2]
labels_full = ["\n".join(np.array(i.split("_"))[[0,3]]) for i in list(focal_data1)]
axs[0].boxplot(focal_data1, notch=True, patch_artist=True, widths=widths,
           boxprops=dict(facecolor=colors[0], color=colors[0]),
           capprops=dict(color=colors[0]),
           whiskerprops=dict(color=colors[0]),
           flierprops=dict(color=colors[0], markeredgecolor=colors[0]),
           medianprops=dict(color=colors[0]))
axs[0].boxplot(focal_data2, positions=np.arange(len(focal_data2.columns))+1.3,
           notch=True, patch_artist=True, widths=widths,
           boxprops=dict(facecolor=colors[1], color=colors[1]),
           capprops=dict(color=colors[1]),
           whiskerprops=dict(color=colors[1]),
           flierprops=dict(color=colors[1], markeredgecolor=colors[1]),
           medianprops=dict(color=colors[1])
           )
axs[0].set_title('High polygenicity')
focal_min = np.concat([focal_data1, focal_data2], axis=1).min(axis=None)
focal_max = np.concat([focal_data1, focal_data2], axis=1).max(axis=None)
focal_range = focal_max-focal_min
axs[0].set_ylim(focal_min - focal_range*0.1,
        focal_max + focal_range*0.2)
# ax.get_xticks()
axs[0].set_xticks(ticks=np.arange(len(focal_data2.columns))+1.125,
           labels=labels_map,
           rotation=0,
           ha='center')
labels2 = axs[0].secondary_xaxis(location=0)
labels2.set_xticks([1.7, 3.7, 5.7, 7.7], labels=labels_model)
labels2.tick_params('x', length=0)
axs[0].set_xlabel(xlabel="\n\nModels")
axs[0].set_ylabel(ylabel="Minimum number of alleles for explaining 80% of\n"+r"positive contribution to local adaptation ( $\it{K_{80,pos}}$)")
trans = axs[0].get_xaxis_transform()
axs[0].plot([0.75,2.5],[-.09,-.09], color="dimgrey", transform=trans, clip_on=False)
axs[0].plot([2.75,4.5],[-.09,-.09], color="dimgrey", transform=trans, clip_on=False)
axs[0].plot([4.75,6.5],[-.09,-.09], color="dimgrey", transform=trans, clip_on=False)
axs[0].plot([6.75,8.5],[-.09,-.09], color="dimgrey", transform=trans, clip_on=False)
axs[0].legend(handles=[patch2, patch1],
           title="",
           loc="upper left")

#Lower plot: low polygenicity
focal_data1 = k80_lowPoly_lowMig
focal_data2 = k80_lowPoly_highMig
labels_map = [i.split("_")[3] for i in list(focal_data1)]
labels_model = ["\n\n"+i.split("_")[0] for i in list(focal_data1)][0::2]
labels_full = ["\n".join(np.array(i.split("_"))[[0,3]]) for i in list(focal_data1)]
axs[1].boxplot(focal_data1, notch=True, patch_artist=True, widths=widths,
           boxprops=dict(facecolor=colors[0], color=colors[0]),
           capprops=dict(color=colors[0]),
           whiskerprops=dict(color=colors[0]),
           flierprops=dict(color=colors[0], markeredgecolor=colors[0]),
           medianprops=dict(color=colors[0]))
axs[1].boxplot(focal_data2, positions=np.arange(len(focal_data2.columns))+1.3,
           notch=True, patch_artist=True, widths=widths,
           boxprops=dict(facecolor=colors[1], color=colors[1]),
           capprops=dict(color=colors[1]),
           whiskerprops=dict(color=colors[1]),
           flierprops=dict(color=colors[1], markeredgecolor=colors[1]),
           medianprops=dict(color=colors[1])
           )
axs[1].set_title('Low polygenicity')
focal_min = np.concat([focal_data1, focal_data2], axis=1).min(axis=None)
focal_max = np.concat([focal_data1, focal_data2], axis=1).max(axis=None)
focal_range = focal_max-focal_min
axs[1].set_ylim(focal_min - focal_range*0.1,
        focal_max + focal_range*0.2)
# ax.get_xticks()
axs[1].set_xticks(ticks=np.arange(len(focal_data2.columns))+1.125,
           labels=labels_map,
           rotation=0,
           ha='center')
labels2 = axs[1].secondary_xaxis(location=0)
labels2.set_xticks([1.7, 3.7, 5.7, 7.7], labels=labels_model)
labels2.tick_params('x', length=0)
axs[1].set_xlabel(xlabel="\n\nModels")
axs[1].set_ylabel(ylabel="Minimum number of alleles for explaining 80% of\n"+r"positive contribution to local adaptation ( $\it{K_{80,pos}}$)")
trans = axs[1].get_xaxis_transform()
axs[1].plot([0.75,2.5],[-.09,-.09], color="dimgrey", transform=trans, clip_on=False)
axs[1].plot([2.75,4.5],[-.09,-.09], color="dimgrey", transform=trans, clip_on=False)
axs[1].plot([4.75,6.5],[-.09,-.09], color="dimgrey", transform=trans, clip_on=False)
axs[1].plot([6.75,8.5],[-.09,-.09], color="dimgrey", transform=trans, clip_on=False)
axs[1].legend(handles=[patch2, patch1],
           title="",
           loc="upper left")
fig.tight_layout()
plt.savefig(figPath+
            "k80_32models_boxplot.png",
            dpi=300)
plt.close()