#!/usr/bin/env bash
#Merge tables 
cd /home/anadem/github/data/tskit_data/output/mutiRuns/k80
touch k80_28model.tab
for file in ./singleModel/*.txt
do 
  paste k80_30model.tab $file > temp.txt
  mv temp.txt k80_28model.tab
done

#Remove the tabs from the begining of each line 
sed -i "s/^[ \t]*//" k80_28model.tab
