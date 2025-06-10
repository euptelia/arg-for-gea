"""
Plot k80 and fst from simulations
tianlin.duan42@gmail.com
2025.03.19
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
#matplotlib.use('Qt5Agg')
import scipy.stats as stats
import matplotlib.patches as mpatches # manually make legends

############################# program #########################################
figPath = "/home/anadem/github/data/tskit_data/figure/multiRuns/"
k80 = pd.read_table("/home/anadem/github/data/tskit_data/output/mutiRuns/k80/k80_32model.tab",
                    sep="\t")

# select columns containing 'lowPoly'
k80_lowPoly = k80.filter(like="lowPoly", axis=1)
k80_highPoly = k80.filter(like="highPoly", axis=1)
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
print(k80_lowPoly.stack().sem())
print(k80_highPoly.stack().sem())

#expansion in highPoly
print(k80.filter(like="a_").mean(axis=None))
print(k80.filter(like="b_").mean(axis=None))
print(k80.filter(like="a_").stack().std())
print(k80.filter(like="b_").stack().std())
print(k80.filter(like="a_").stack().sem())
print(k80.filter(like="b_").stack().sem())
# print(k80_highPoly.filter(like="a_").mean(axis=None))
# print(k80_highPoly.filter(like="b_").mean(axis=None))
# print(k80_highPoly.filter(like="a_").stack().sem())
# print(k80_highPoly.filter(like="b_").stack().sem())


#Migration
print(k80_lowMig.stack().mean(axis=None))
print(k80_highMig.stack().mean(axis=None))
print(k80_lowMig.stack().sem()) #
print(k80_highMig.stack().sem()) #
# print(k80_highPoly_lowMig.stack().mean(axis=None))
# print(k80_highPoly_highMig.stack().mean(axis=None))
# print(k80_highPoly_lowMig.stack().sem()) #
# print(k80_highPoly_highMig.stack().sem()) #
# print(k80_lowPoly_lowMig.stack().mean(axis=None))
# print(k80_lowPoly_highMig.stack().mean(axis=None))
# print(k80_lowPoly_lowMig.stack().sem()) #
# print(k80_lowPoly_highMig.stack().sem()) #

#cline/patchy
print(k80.filter(like="Cline").mean(axis=None)) #
print(k80.filter(like="Patchy").mean(axis=None)) #
print(k80.filter(like="Cline").stack().sem()) #
print(k80.filter(like="Patchy").stack().sem()) #

# print(k80_highPoly.filter(like="Cline").mean(axis=None)) #
# print(k80_highPoly.filter(like="Patchy").mean(axis=None)) #
# print(k80_highPoly.filter(like="Cline").stack().sem()) #
# print(k80_highPoly.filter(like="Patchy").stack().sem()) #


#One way ANOVA
f_oneway(k80.filter(like="a_").stack(),
         k80.filter(like="b_").stack())
f_oneway(k80_lowMig.stack(),
         k80_highMig.stack())
f_oneway(k80.filter(like="Cline").stack(),
         k80.filter(like="Patchy").stack())

#Only high poly
f_oneway(k80_highPoly.filter(like="a_").stack(),
         k80_highPoly.filter(like="b_").stack())
f_oneway(k80_highPoly_lowMig.stack(),
         k80_highPoly_highMig.stack())
f_oneway(k80_highPoly.filter(like="Cline").stack(),
         k80_highPoly.filter(like="Patchy").stack())



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
labels1 = ["\n".join(np.array(i.split("_"))[[0,3]]) for i in list(focal_data1)]
labels2 = ["\n".join(np.array(i.split("_"))[[0,3]]) for i in list(focal_data2)]
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
plt.ylim(focal_min - focal_range*0.1,
        focal_max + focal_range*0.3)
# ax.get_xticks()
axs[0].set_xticks(ticks=np.arange(len(focal_data2.columns))+1.125,
           labels=labels2,
           rotation=20,
           ha='center')
axs[0].set_xlabel(xlabel="Models")
axs[0].set_ylabel(ylabel="Minimum number of alleles for explaining 80% of\n"+r"positive contribution to local adaptation ( $\it{K_{80}(positive)}$)")
axs[0].legend(handles=[patch2, patch1],
           title="",
           loc="upper left")

#Lower plot: low polygenicity
focal_data1 = k80_lowPoly_lowMig
focal_data2 = k80_lowPoly_highMig
labels1 = ["\n".join(np.array(i.split("_"))[[0,3]]) for i in list(focal_data1)]
labels2 = ["\n".join(np.array(i.split("_"))[[0,3]]) for i in list(focal_data2)]
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
plt.ylim(focal_min - focal_range*0.1,
        focal_max + focal_range*0.3)
# ax.get_xticks()
axs[1].set_xticks(ticks=np.arange(len(focal_data2.columns))+1.125,
           labels=labels2,
           rotation=20,
           ha='center')
axs[1].set_xlabel(xlabel="Models")
axs[1].set_ylabel(ylabel="Minimum number of alleles for explaining 80% of\n"+r"positive contribution to local adaptation ( $\it{K_{80}(positive)}$)")
axs[1].legend(handles=[patch2, patch1],
           title="",
           loc="upper left")
fig.tight_layout()
plt.savefig(figPath+
            "k80_32models_boxplot.png",
            dpi=300)
plt.close()