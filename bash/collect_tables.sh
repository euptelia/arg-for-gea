#!/usr/bin/env bash


#Move files together
find /home/anadem/github/data/tskit_data/output/mutiRuns -type f -mindepth 2 -exec mv -i -- {} /home/anadem/github/data/tskit_data/output/k80 \;

#Replace long names with shorter ones
rename 's/Continuous_nonWF_M3b_glacialHistoryOptimum0/M3b/' /home/anadem/github/data/tskit_data/output/k80/renamed/*.txt
rename 's/Continuous_nonWF_M3a_glacialHistoryOptimum0/M3a/' /home/anadem/github/data/tskit_data/output/k80/renamed/*.txt
rename 's/Continuous_nonWF_M2b_glacialHistoryOptimum0/M2b/' /home/anadem/github/data/tskit_data/output/k80/renamed/*.txt
rename 's/Continuous_nonWF_M2a_glacialHistoryOptimum0/M2a/' /home/anadem/github/data/tskit_data/output/k80/renamed/*.txt

#polygeniciy
rename 's/mu1.0e-10_sigmaM0.1/lowPoly/' /home/anadem/github/data/tskit_data/output/k80/renamed/*.txt
rename 's/mu1.0e-08_sigmaM0.01/highPoly/' /home/anadem/github/data/tskit_data/output/k80/renamed/*.txt

#migration and the rest
rename 's/sigmaW0.4_sigmaD0.03_mateD0.12_K17000_r1.0e-07_tick110000_200runs/lowMig/' /home/anadem/github/data/tskit_data/output/k80/renamed/*.txt
rename 's/sigmaW0.4_sigmaD0.06_mateD0.15_K17000_r1.0e-07_tick110000_200runs/highMig/' /home/anadem/github/data/tskit_data/output/k80/renamed/*.txt


#Replace the first line with file names
cd /home/anadem/github/data/tskit_data/output/k80/renamed/
for file in ./*.txt
do 
  header=$(basename $file .txt)
  sed -i "1 s/.*/${header}/" $file
done

#Merge tables 
touch k80_22model.tab
for file in ./*.txt
do 
  paste k80_22model.tab $file > temp.txt
  mv temp.txt k80_22model.tab
done

#Remove the tabs from the begining of each line 
sed -i "s/^[ \t]*//" k80_22model.tab
