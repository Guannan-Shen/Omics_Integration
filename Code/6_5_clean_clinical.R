## clean the clinical Data 
## with and a list of samples 
'%nin%' <- Negate('%in%')
options(stringsAsFactors = F)

library(readxl)
library(tidyverse)
library(magrittr)
setwd("~/Documents/gitlab/Omics_Integration/DataRaw/hiv_infected_un/")
getwd()


# cli2 <- read.xlsx("clin_microbiome.xlsx")

## subset samples (ID Lib) and rescale 
rescaled_cli <- function(){
  ############## 26 samples ###########
  small_lib = c("MIHIV998")
  ## shared sample size across datasets
  shared_sam = c( "MIHIV124", "MIHIV132", "MIHIV138", "MIHIV154", "MIHIV178", "MIHIV255", 
                  "MIHIV278", "MIHIV286", "MIHIV323", "MIHIV361", "MIHIV391", "MIHIV404", 
                  "MIHIV428", "MIHIV493", "MIHIV582", "MIHIV594", "MIHIV648", "MIHIV683", 
                  "MIHIV708", "MIHIV716", "MIHIV819", "MIHIV825", "MIHIV839", "MIHIV914", 
                  "MIHIV947", "MIHIV972", "MIHIV998")
  final_sam = shared_sam[ shared_sam %nin% small_lib]
  
  # get clinical data
  cli1 = read.xlsx("clinical_ready.xlsx")
  # subset
  data = cli1 %>% dplyr::filter(cli1$pid %in% final_sam)
  # rescale to mean 0 and variance 1, from age to CD4Tcells % viable CD45+
  df = cbind(data[, 1:4], apply(data[, 5:25], 2, scale) )
  # make the clinical data has the same order of samples 
  # genelist contains Gene_ID, Symbol
  # clean and transform transcriptome data, subset of genes
  source("~/Documents/gitlab/Omics_Integration/Code/6_5_clean_transcriptome.R")
  # rlog transform or not
  rnaseq = rescaled_rna(rlog = T)[[1]] %>% rownames_to_column("pid") %>% select(pid)
  # plyr::join keep the order of the left dataset
  results = plyr::join(rnaseq, df, by = "pid" )
  colnames(results)[1] = "ID"
  
  return(results)
}

