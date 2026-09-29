#!/usr/bin/env bash
# intermediate mutational polygenicity
#for i in {1..50} 
#do
#    printf "Start the %s/50 runs of M3a, cline, intermediate Poly, high Mig\n" $i
#    slim /home/anadem/github/arg-for-gea/slim/second_revision/intermediate_polygenicity/contiuous_nonWF_M3a_glacialHistory_clineMap_recurrentChange.slim
#done

## intermediate mutational polygenicity
#for i in {1..50} 
#do
#    printf "Start the %s/50 runs of M3b, cline, intermediate Poly, high Mig\n" $i
#    slim /home/anadem/github/arg-for-gea/slim/second_revision/intermediate_polygenicity/contiuous_nonWF_M3b_glacialHistory_clineMap_recurrentChange.slim
#done

# intermediate mutational polygenicity
#for i in {1..50} 
#do
#    printf "Start the %s/50 runs of M3a, cline, intermediate Poly, high Mig\n" $i
#    slim /home/anadem/github/arg-for-gea/slim/second_revision/intermediate_polygenicity/contiuous_nonWF_M3a_glacialHistory_clineMap_recurrentChange_intermediate2.slim
#done

# intermediate mutational polygenicity
for i in {1..200} 
do
    printf "Start the %s/200 runs of M3b, cline, intermediate Poly2, high Mig\n" $i
    slim /home/anadem/github/arg-for-gea/slim/second_revision/intermediate_polygenicity/contiuous_nonWF_M3b_glacialHistory_clineMap_recurrentChange_intermediate2.slim
done
