#!/usr/bin/env bash
# Plot LF_pop for one replicate per model
inPath=/media/anadem/PortableSSD/arg4gea_data/slim_data
for folder in ${inPath}/M*/tick110000/
do
    file1=$(ls $folder | head -n 1)
    echo $file1
    python3 /home/anadem/github/arg-for-gea/python/second_revision/lf_by_pop.py -i ${inPath}/M*/tick110000/${file1}
done
