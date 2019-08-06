rm(list = ls())
library(openxlsx)
library(tidyverse)
library(magrittr)
library(SmCCNet)

options(stringsAsFactors = F)
options(dplyr.width = Inf)
setwd("~/Documents/gitlab/Omics_Integration/")
getwd()

# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
# small functions, %nin%, %||%, isnothing
source( paste0(dir, "Code/small_fun.R") )

########## import dataset ################
# source to get the load_filtered_micro_level function to get clr of RA
source( paste0(dir, "Code/5_29_Generate_filtered_Data_Microbiome.R") )
# clean and transform transcriptome data, subset of genes
source( paste0(dir, "Code/6_5_clean_transcriptome.R") )
# clinical data 
source( paste0(dir, "Code/6_5_clean_clinical.R") )
# diagnostic plots and tables
source( paste0(dir, "Code/ref_plots.R") )

###### load in ####
#### phenotype contains ID
clin <- rescaled_cli() 
# get hints ideas
# the final clinical parameter to use
# where na in LPS
n_na <-  which(is.na( clin$LPS ))
LPS <- clin %>% select(LPS) %>% na.omit() 

# import gene lists #
isgs <- as.data.frame(read.delim( paste0(dir, "DataRaw/hiv_infected_un/coreISG")) )
rna_isgs <- rescaled_rna(isgs, rlog = T)
# rlog transform or not
isgs_rlog <- rna_isgs[[1]]
# names
colnames(isgs_rlog) <- rna_isgs[[2]]$Symbol

#### Microbiome
micro_data <- load_filtered_micro_level_samples("genus",  prevalence = 40, RA = 2, wd = "Ubuntu")
micro_clr <- micro_data[[2]] %>% as.data.frame()

# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)


