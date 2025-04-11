#!/usr/bin/env bash
basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection/Continuous_nonWF_M2b_glacialHistoryOptimum0_*/timeSeries/"
for folder in $(echo $basePath)
do
  echo $folder
  python3 /home/anadem/github/arg-for-gea/python/timeSeries_multiTimes_multiRuns_crutches.py -i $folder
done

