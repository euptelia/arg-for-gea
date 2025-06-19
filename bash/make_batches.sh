#!/usr/bin/env bash
nBatch=4
# inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M0b_highPoly_highMig/"
# inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M2b_highPoly_lowMig_patchyMap/"
# inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/historical_optimum_0/M2b_highPoly_lowMig_clineMap/"
# inPath="/home/anadem/github/data/slim_data/glacial_history/20241127variedrefugia/m2b_patchy_highMig_variedRefugia/tick102000/"
inPath="/home/anadem/github/data/slim_data/glacial_history/recurrentChange/M3a_smallLowVm_highMig_patchyMap/timeSeries/"
cd ${inPath}
mkdir other_times
mv *tick100000* other_times
mv *tick101200* other_times
mv *tick101400* other_times
mv *tick101600* other_times
mv *tick101800* other_times
nfiles=$(ls *.trees | wc -l)
nPerbatch=$((nfiles / nBatch))
for (( i=1; i<=${nBatch}; i++ ))
do
  mkdir batch${i}
  mv `ls *.trees | head -n ${nPerbatch}` ./batch${i}/
done

