#!/usr/bin/env bash
#basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection"
#basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection/*M2a*mu1.0e-10_sigmaM0.1_sigmaW0.4_sigmaD0.03_mateD0.12_K17000_r1.0e-07*"
#for folder in $(ls $basePath)
#do
#  echo $basePath/$folder/tick110000/
#  python3 /home/anadem/github/arg-for-gea/python/alleleAge_multiRuns_crutches.py -i $basePath/$folder/tick110000/
#done

#!/usr/bin/env bash
#basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection/Continuous_nonWF_*"
#sleep 5h
basePath="/home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection/*M2a*mu1.0e-10_sigmaM0.1_sigmaW0.4_sigmaD0.03_mateD0.12_K17000_r1.0e-07*"
for folder in $basePath
do
  echo $folder/tick110000/
  python3 /home/anadem/github/arg-for-gea/python/alleleAge_multiRuns_crutches.py -i $folder/tick110000/
done

#python3 /home/anadem/github/arg-for-gea/python/alleleAge_multiRuns_crutches.py -i /home/anadem/github/data/tskit_data/output/table/realistic_fpr_comparisons/selection/#Continuous_nonWF_M3b_glacialHistoryOptimum0_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.03_mateD0.12_K17000_r1.0e-07/tick110000/


