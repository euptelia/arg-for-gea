#!/usr/bin/env bash
# basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection/Continuous_nonWF_*/timeSeries/"
# basePath="/media/anadem/PortableSSD/arg4gea_data/tskit_data/output/table/Continuous_nonWF_*/timeSeries/"
# basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection/Continuous_nonWF_M3b_glacialHistoryOptimum0_clineMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15*/timeSeries/"
basePath="/media/anadem/golden5tb/arg4gea_data/tskit_data/timeSeries/Continuous_nonWF_*/timeSeries/"
for folder in $basePath
do
  echo $folder
  #mkdir $folder/tick100000
  #mv $folder/*_tick100000_* $folder/tick100000
  python3 /home/anadem/github/arg-for-gea/python/timeSeries_multiTimes_multiRuns_crutches.py -i $folder
done

