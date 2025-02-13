#!/usr/bin/env bash
basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection"
for folder in $(ls $basePath)
do
  echo $basePath/$folder/tick110000/
  python3 /home/anadem/github/arg-for-gea/python/alleleAge_multiRuns_crutches.py -i $basePath/$folder/tick110000/
done

