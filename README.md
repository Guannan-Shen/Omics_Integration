# Omics Integration
## Background
Explore state of the art methods in Omics data and Phenotype integration, especially the integration and network analysis of **RNAseq, Microbiome and HIV-related clinical parameters**.  

Since I have a statistics background, the marjority of this analysis will be done in R.  

Typically, these methods could be very slow to run. For example, SmCCNet could take 1-2 hours to run via CPU.

## Methods
1. [SmCCNet](https://cran.rstudio.com/web/packages/SmCCNet/index.html): Sparse Multiple Canonical Correlation Network Analysis Tool.  
2. [mixOmics](https://www.bioconductor.org/packages/release/bioc/html/mixOmics.html): Omics Data Integration Project.  
3. [r.jive](https://cran.r-project.org/web/packages/r.jive/r.jive.pdf): Performs the JIVE decomposition on a list of data sets when the data share a dimension, returning low-rank matrices that capture the joint and individual structure of the data.  

## Folders, files and datasets
1. All datasets are saved locally under DataRaw/hiv_infected_un/. The raw datasets were loaded in using 8_5_testing_dataset.R in Code/.  
2. Workflow using scripts after date 8/5/2019 (8_5). 
3. The DataProcessed directory contains all results and figures generated for different cases using the pipeline. The main pipeline scripts are get_l1l2.R, run_SmCCnet.R and after_running.R. 
