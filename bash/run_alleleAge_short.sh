#!/usr/bin/env bash
batchNumber=1
inPath="/home/tianlin/ubc/data/slim/20241127variedrefugia/m2b_patchy_highMig_variedRefugia/tick110000/batch${batchNumber}"
j=1
echo $inPath
for file in ${inPath}/*.trees
do
    printf "file %s / 50 \n" $j 
    printf "Start making a table with \n %s" $file
    python3 /home/tianlin/ubc/arg-for-gea/python/alleleAge_slimHistory_crutches.py -i ${file} -p 0
    j=$((${j}+1))
done

