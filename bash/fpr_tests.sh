#!/usr/bin/env bash
#Second revision: FPR tests
# Only neutral chromosome
basePath="/media/anadem/PortableSSD/arg4gea_data/tskit_data/tick110000/Continuous_nonWF*"
for folder in $basePath
do
  echo $folder/tick110000/
  python3 /home/anadem/github/arg-for-gea/python/second_revision/calculate_fpr.py -i $folder/tick110000/
done

basePath="/media/anadem/PortableSSD/arg4gea_data/tskit_data/wrong_fitness/neutral/Continuous_nonWF*"
for folder in $basePath
do
  echo $folder/tick110000/
  python3 /home/anadem/github/arg-for-gea/python/second_revision/calculate_fpr_m0.py -i $folder/tick110000/
done
