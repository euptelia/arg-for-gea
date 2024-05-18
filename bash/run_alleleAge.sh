#!/usr/bin/env bash
#for file in /home/tianlin/Documents/github/data/slim_data/m2b_mu1e-7_m0.1_w0.4/test/*.trees
#do
#    modelName=$(basename $file .trees)
#    python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_neutralMut.py -t $modelName
#done

#batchNumber=1
#inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/highPoly_highMig_clineMap/batch${batchNumber}"
#for file in ${inPath}/*.trees
#do
#    modelName=$(basename $file .trees)
#    python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -t batch${batchNumber}/$modelName
#done

# batchNumber=4
# inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_highPoly_highMig_patchyMap/batch${batchNumber}"
# for file in ${inPath}/*.trees
# do
#     python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -i ${file}
# done

# batchNumber=4
# inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M0b_highPoly_highMig/batch${batchNumber}"
# for file in ${inPath}/*.trees
# do
#     python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -i ${file}
# done

# batchNumber=1
# inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_highPoly_lowMig_patchyMap/batch${batchNumber}"
# for file in ${inPath}/*.trees
# do
#     python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -i ${file}
# done

#batchNumber=1
#inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_highPoly_lowMig_clineMap/batch${batchNumber}"
#for file in ${inPath}/*.trees
#do
#    python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -i ${file}
#done

#batchNumber=1
#inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_lowPoly_lowMig_clineMap/batch${batchNumber}"
#for file in ${inPath}/*.trees
#do
#    python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -i ${file}
#done

#batchNumber=1
#inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_lowPoly_lowMig_patchyMap/batch${batchNumber}"
#for file in ${inPath}/*.trees
#do
#    python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -i ${file}
#done

# batchNumber=1
# inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_highPoly_lowMig_patchyMap/batch${batchNumber}"
# for file in ${inPath}/*.trees
# do
#     python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -i ${file}
# done

# batchNumber=1
# inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_lowPoly_highMig_patchyMap/batch${batchNumber}"
# for file in ${inPath}/*.trees
# do
#    python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory_temp2.py -i ${file}
# done

# batchNumber=6
# inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_highPoly_highMig_clineMap/batch${batchNumber}"
# for file in ${inPath}/*.trees
# do
#    python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory_temp2.py -i ${file}
# done

#batchNumber=1
#inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_lowPoly_highMig_clineMap/batch${batchNumber}"
#for file in ${inPath}/*.trees
#do
#   python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory_temp2.py -i ${file}
#done

#batchNumber=4
#inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/historical_optimum_0/M2b_highPoly_highMig_clineMap/batch${batchNumber}"
#for file in ${inPath}/*.trees
#do
#   python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory_temp2.py -i ${file}
#done

#batchNumber=4
#inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/historical_optimum_0/M2b_highPoly_lowMig_clineMap/batch${batchNumber}"
#for file in ${inPath}/*.trees
#do
#   python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -i ${file}
#done

# batchNumber=4
# inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/historical_optimum_0/M2b_lowPoly_highMig_clineMap/batch${batchNumber}"
# for file in ${inPath}/*.trees
# do
#    python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -i ${file}
# done

batchNumber=4
inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/historical_optimum_0/M2b_lowPoly_lowMig_clineMap/batch${batchNumber}"
for file in ${inPath}/*.trees
do
   python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_slimHistory.py -i ${file}
done

