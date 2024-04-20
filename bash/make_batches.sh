#!/usr/bin/env bash
nBatch=4
inPath="/home/tianlin/Documents/github/data/slim_data/glacial_history/M0b_highPoly_highMig/"
cd ${inPath}
nfiles=$(ls *.trees | wc -l)
nPerbatch=$((nfiles / nBatch))
for (( i=1; i<=${nBatch}; i++ ))
do
  mkdir batch${i}
  mv `ls *.trees | head -n ${nPerbatch}` ./batch${i}/
done

