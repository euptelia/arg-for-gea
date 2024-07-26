"""
For crushed jobs, check which input files have been used and move them into a
folder named "used".
"""
############################# modules #########################################
import os # mkdir
# ############################# options #############################
# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('-i', '--inPath',
#                     help='Main path to the input files',
#                     type=str)
# parser.add_argument('-o', '--outPath',
#                     help='Path to the output files',
#                     type=int, default=1)
# args = parser.parse_args()
#
# ############################# program #########################################
# inPath = args.inPath
# outPath = args.outPath
# inPath = "/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_highPoly_lowMig_patchyMap/"
# outPath = "/home/tianlin/Documents/github/data/tskit_data/output/table/Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.03_mateD0.12/"
# inPath = "/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_lowPoly_lowMig_patchyMap/"
# outPath = "/home/tianlin/Documents/github/data/tskit_data/output/table/Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-09_sigmaM0.1_sigmaW0.4_sigmaD0.03_mateD0.12/"
# inPath = "/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_lowPoly_lowMig_clineMap/"
# outPath = "/home/tianlin/Documents/github/data/tskit_data/output/table/Continuous_nonWF_M2b_glacialHistory_clineMap_mu1.0e-09_sigmaM0.1_sigmaW0.4_sigmaD0.03_mateD0.12/"
# inPath = "/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_highPoly_lowMig_clineMap/"
# outPath = "/home/tianlin/Documents/github/data/tskit_data/output/table/Continuous_nonWF_M2b_glacialHistory_clineMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.03_mateD0.12/"
# inPath = "/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_lowPoly_highMig_patchyMap/"
# outPath = "/home/tianlin/Documents/github/data/tskit_data/output/table/Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-09_sigmaM0.1_sigmaW0.4_sigmaD0.06_mateD0.15/"
inPath = "/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_lowPoly_highMig_clineMap/"
outPath = "/home/tianlin/Documents/github/data/tskit_data/output/table/Continuous_nonWF_M2b_glacialHistory_clineMap_mu1.0e-09_sigmaM0.1_sigmaW0.4_sigmaD0.06_mateD0.15/tick101000/"


usedNames = ["_".join(f.split("_")[0:12]) #Check this before use!
             for f in os.listdir(outPath)
             if os.path.isfile(os.path.join(outPath, f))]

inFolders = [os.path.join(inPath, f) for f in os.listdir(inPath)
             if os.path.isdir(os.path.join(inPath, f))]
useFilePath = os.path.join(inPath, "used")
if not os.path.exists(useFilePath):
    os.makedirs(useFilePath)
for batch in inFolders:
    for f in os.listdir(batch):
        if f.split(".tree")[0] in usedNames:
            # print(os.path.join(batch, f))
            # print(os.path.join(useFilePath, f))
            os.rename(os.path.join(batch, f), os.path.join(useFilePath, f))
