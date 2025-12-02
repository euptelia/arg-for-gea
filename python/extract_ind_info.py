"""
For all .trees files in a folder,
Extract:
    1. Number of individuals
    2. Individual ages
Save results in tables
About 5 min for
tianlin.duan42@gmail.com
Last modified
2025.11.10
"""
############################# modules #########################################
import tskit
import os # mkdir
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
inPath = args.inPath
outPath = args.outPath
fileList = glob.glob(inPath + "*.trees")
model_name = fileList[0].split("/")[-1].split(".trees")[0]

#Translate long model names to short model names
if "sigmaD0.06_mateD0.15" in model_name:
    migName = "HighMig"
elif "sigmaD0.03_mateD0.12" in model_name:
    migName = "LowMig"
else:
    migName = ""
if "_clineMap_" in model_name:
    mapName = "Cline"
elif "_patchyMap_" in model_name:
    mapName = "Patchy"
else:
    mapName = ""
if "_sigmaM0.01_" in model_name:
    mutName = "HighPoly"
elif "_sigmaM0.1_" in model_name:
    mutName = "LowPoly"
else:
    mutName = ""
demoName_ori = model_name.split("_")[2]
name_change = {"M2a":"MSelCon", "M2b":"MSelExp", "M3a":"MRecCon", "M3b":"MRecExp"}
demoName = name_change[demoName_ori]
shortName = "_".join([demoName, mutName, migName, mapName])

#Create directories
for folder in ["/ind_age/", "/n_ind/"]:
    if not os.path.exists(outPath+folder):
        os.makedirs(outPath+folder)
output_indAge = outPath + "/ind_age/" + shortName +"_indAge.tab"
output_n_ind = outPath + "/n_ind/" + shortName +"_n_ind.tab"
fout_indAge = open(output_indAge, "w")
fout_n_ind = open(output_n_ind, "w")
fout_indAge.write(shortName+"\n")
fout_n_ind.write(shortName+"\n")

#Loading individual .tree files
for file in fileList:
    # print("Fst calculation started")
    # print(time.ctime())
    file_name = file.split("/")[-1]
    # t1 = time.time()
    # Tree-sequence file from SLiM
    ts = tskit.load(file)

    # Number of diploid individuals
    N = ts.num_individuals
    fout_n_ind.write(str(N) + "\n")
    # Age of individuals
    for ind in ts.individuals():
        fout_indAge.write(str(ind.metadata['age'])+"\n")
fout_indAge.close()
fout_n_ind.close()
