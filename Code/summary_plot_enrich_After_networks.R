rm(list = ls())
library(openxlsx)

library(magrittr)
library(SmCCNet)
library(tidyverse)
options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()

# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
# small functions, %nin%, %||%, isnothing
source( paste0(dir, "Code/small_fun.R") )
# ref plots and correlation 
source( paste0(dir, "Code/ref_plots.R") )
# get module 0 and PCA
source( paste0(dir, "Code/corr_pheno.R") )
# Optimize edge cut 
source( paste0(dir, "Code/edge_cut.R") )
source( paste0(dir, "Code/put_together.R") )
source( paste0(dir, "Code/enrichment.R") )

source( paste0(dir, "Code/pca_getPC1.R") )

# load all datasets 
source( paste0(dir, "Code/8_5_testing_dataset.R") )

####################

################# load abar, modules, Ws, the product of SmCCNet ###############
CVDir <- "CD14_Outlier3_Global_100_50_Genus_3_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))
dim(abar)

