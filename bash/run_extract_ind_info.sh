#!/usr/bin/env bash
#Run extract_ind_info.py for multiple folders
basePath="/media/anadem/PortableSSD/arg4gea_data/slim_data/*/tick110000/"
# basePath="/media/anadem/PortableSSD/arg4gea_data/slim_data/M2b_smallLowVm*/tick110000/"
for folder in $basePath
do
  echo $folder
  cd $folder
  python3 /home/anadem/github/arg-for-gea/python/extract_ind_info.py -i $folder -o /home/anadem/github/data/tskit_data/stats/test
done
