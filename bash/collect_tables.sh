#!/usr/bin/env bash
#Merge k80 tables 
#cd /home/anadem/github/data/tskit_data/output/mutiRuns/k80
cd /home/anadem/github/data/tskit_data/stats/k80
touch k80_32model.tab
for file in ./singleModel/*.txt
do 
  paste k80_32model.tab $file > temp.txt
  mv temp.txt k80_32model.tab
done

#Remove the tabs from the begining of each line 
sed -i "s/^[ \t]*//" k80_32model.tab

#Move tick110000 figures together
cd /home/anadem/github/data/tskit_data/figure/multiRuns   
mkdir freqColor ageColor k80 lf fpr
mv /home/anadem/github/data/tskit_data/figure/multiRuns/*/*sizeColor.png sizeColor
mv /home/anadem/github/data/tskit_data/figure/multiRuns/*/*_LF_by_age_bin_positve_lfRelative2EachRun_freqColor.png freqColor
mv /home/anadem/github/data/tskit_data/figure/multiRuns/*/*_LF_by_size_bin_positve_lfRelative2PositiveLFEachRun_freqColor.png lfBySize_freqColor
mv /home/anadem/github/data/tskit_data/figure/multiRuns/*/*GEAp1e-10_maf0.05_equalInterval.png fpr       

#Collect number of individuals in each simulation
cd /home/anadem/github/data/tskit_data/stats/n_ind
touch n_ind_32model.tab
for file in ./M*_n_ind.tab
do 
  paste n_ind_32model.tab $file > temp.txt
  mv temp.txt n_ind_32model.tab
done
#Remove the tabs from the begining of each line 
sed -i "s/^[ \t]*//" n_ind_32model.tab

#Collect the summary of fst and lf
cd /home/anadem/github/data/tskit_data/stats/summary
head -n 1 MRecCon_HighPoly_HighMig_Cline_averageFst_averageLF_averageAge.tab > fst_lf_32model.tab
for file in ./M*averageFst_averageLF*.tab
do 
  tail -n +2 $file >> fst_lf_32model.tab
done

#Gather figures together
cd /home/anadem/github/data/tskit_data/figure/multiTimes/multiRun
mkdir k80 lf mutTime freq size
mv /home/anadem/github/data/tskit_data/figure/multiTimes/multiRun/*/*_lfByMutTime_20Bin_stack.png mutTime
mv /home/anadem/github/data/tskit_data/figure/multiTimes/multiRun/*/*_20RunsMean_lfByFreq10Bin_stack.png freq
mv /home/anadem/github/data/tskit_data/figure/multiTimes/multiRun/*/*_20RunsMean_lfBysize20Bin_stack.png size
mv /home/anadem/github/data/tskit_data/figure/multiTimes/multiRun/*/*_20RunsMedian_allelicLF.png lf
mv /home/anadem/github/data/tskit_data/figure/multiTimes/multiRun/*/*_20RunsMedian_k80net.png k80

cd /home/anadem/github/data/tskit_data/figure/multiTimes/multiRun
mkdir k80 lf mutTime freq size
