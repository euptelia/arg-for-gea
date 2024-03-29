#!/usr/bin/env bash
for file in /home/tianlin/Documents/github/data/slim_data/m2b_mu1e-7_m0.1_w0.4/test/*.trees
do
    modelName=$(basename $file .trees)
    python3 /home/tianlin/Documents/github/arg-for-gea/python/alleleAge_neutralMut.py -t $modelName
done
