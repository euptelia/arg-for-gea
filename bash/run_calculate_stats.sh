#!/usr/bin/env bash
#Run calculate_stats.py for multiple folders
#basePath="/media/anadem/PortableSSD/arg4gea_data/slim_data/M2a_largeLowVm_highMig_clineMap/tick110000/"
#for folder in $(echo $basePath)
#do
#  echo $folder
#  cd $folder
#  #Move .trees files together
#  for subfolder in ${folder}*/
#  do
#    mv $subfolder*trees .
#  done
#  python3 /home/anadem/github/arg-for-gea/python/calculate_stats.py -i $folder -o /home/anadem/github/data/tskit_data/stats
#done

basePath="/media/anadem/PortableSSD/arg4gea_data/slim_data/*/tick110000/"
for folder in $basePath
do
  echo $folder
  cd $folder
  python3 /home/anadem/github/arg-for-gea/python/calculate_stats.py -i $folder -o /home/anadem/github/data/tskit_data/stats
done
