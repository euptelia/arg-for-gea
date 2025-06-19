#!/usr/bin/env bash
#basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection/Continuous_nonWF_*/timeSeries/"
#basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection/Continuous_nonWF_M2a_*/timeSeries/"
#basePath="/media/anadem/PortableSSD/arg4gea_data/tskit_data/Continuous_nonWF_M3a_glacialHistoryOptimum0_patchyMap*/timeSeries/"
basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/Continuous_nonWF_M3a_glacialHistoryOptimum0_patchyMap_mu1.0e-10_*/timeSeries/"
for folder in $(echo $basePath)
do
  echo $folder
  python3 /home/anadem/github/arg-for-gea/python/timeSeries_multiTimes_multiRuns_crutches.py -i $folder
done

