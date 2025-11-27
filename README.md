# arg-for-gea  
SLiM simulations and Python scripts for exploring the genetic architecture of local adaptation and its implications for genotype-environment association (GEA) analysis, with a focus on the effects of allele age. 
Tree sequences are recorded during SLiM simulations and processed by Python scripts.  

*Note: Models have been renamed in reports:  
    Scripts            Reports  
    M0a  --------->    MNeuCon  
    M0b  --------->    MNeuExp  
    M2a  --------->    MSelCon  
    M2a  --------->    MSelExp  
    M3a  --------->    MRecCon  
    M3b  --------->    MRecExp  

### Single-sample simulations  
arg-for-gea/slim/glacial_history/  

#### No selection (M0)  
    contiuous_nonWF_M*.slim  -> alleleAge_slimHistory_m0.py -> alleleAge_multiRuns_m0.py  
#### With selection (M2,M3)  
    contiuous_nonWF_M*.slim  -> alleleAge_slimHistory.py -> alleleAge_multiRuns.py  

### Time-series-sample simulations 
arg-for-gea/slim/glacial_history/time_series/  

    contiuous_nonWF_M*_timeSeries\_*.slim -> alleleAge_slimHistory_noNeuMut.py -> timeSeries_multiTimes_multiRuns.py  

### Calculating statistics
    arg-for-gea/slim/glacial_history/contiuous_nonWF_M*.slim -> calculate_stats.py
