# arg-for-gea  
SLiM simulations and Python scripts for exploring the genetic architecture of local adaptation and its implications for genotype-environment association (GEA) analysis, with a focus on the effects of allele age. 
Tree sequences are recorded during SLiM simulations and processed by Python scripts.  

*Note: Models have been renamed in reports:  
    Scripts  --------->    Reports  
    M0a  --------->    MNeuCon  
    M0b  --------->    MNeuExp  
    M2a  --------->    MSelCon  
    M2a  --------->    MSelExp  
    M3a  --------->    MRecCon  
    M3b  --------->    MRecExp  

### Single-sample simulations  
arg-for-gea/slim/glacial_history/  

#### No selection (M0)  
![\*.slim  -> alleleAge_slimHistory_m0.py -> alleleAge_multiRuns_m0.py](docs/scripts_flow_chart3_1.1.png)  
      
#### With selection (M2,M3)  
![\*.slim  -> alleleAge_slimHistory.py -> alleleAge_multiRuns.py](docs/scripts_flow_chart1_1.1.png)  
 

### Time-series-sample simulations 
arg-for-gea/slim/glacial_history/time_series/  
![\*.slim -> alleleAge_slimHistory_noNeuMut.py -> timeSeries_multiTimes_multiRuns.py](ocs/scripts_flow_chart2_1.1.png)  
