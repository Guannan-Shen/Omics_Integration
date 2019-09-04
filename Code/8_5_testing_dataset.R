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
LPS <- clin %>% select(LPS) # %>% na.omit() 

# check other clinical parameters
sum(which(is.na( clin$CD14 )) )
which(is.na( clin$`IL-6`))
which(is.na( clin$CRP))
sum(which(is.na( LPS )) )
CD14 <- clin %>% select(CD14)

# import gene lists #
isgs <- as.data.frame(read.delim( paste0(dir, "DataRaw/hiv_infected_un/coreISG")) )
rna_isgs <- rescaled_rna(isgs, rlog = T)
# rlog transform or not
isgs_rlog <- rna_isgs[[1]]
# names
colnames(isgs_rlog) <- rna_isgs[[2]]$Symbol

# ifn-beta
genesbeta <- as.data.frame(read.delim( paste0(dir, "DataRaw/hiv_infected_un/genesbeta" )) )
rna_genesbeta <- rescaled_rna(genesbeta, rlog = T)
# rlog transform or not
genesbeta_rlog <- rna_genesbeta[[1]]
# names
colnames(genesbeta_rlog) <- rna_genesbeta[[2]]$Symbol

#### Microbiome
micro_data <- load_filtered_micro_level_samples("genus",  prevalence = 40, RA = 2, wd = "Ubuntu")
micro_clr <- micro_data[[2]] %>% as.data.frame()

# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)

### hints of weights genesbeta ##########
# diagnostic plots and tables
source( paste0(dir, "Code/ref_plots.R") )
df1 = get_corr(mibi[-n_na, ], genesbeta_rlog[-n_na, ], "Microbiome:Beta-ISGs")
# stats::cor(clin$LPS, mibi, use = "pairwise.complete.obs", method = "pearson")
df2 = get_corr(LPS[-n_na], mibi[-n_na, ], "LPS:Microbiome")
df3 = get_corr(LPS[-n_na], genesbeta_rlog[-n_na, ], "LPS:Beta-ISGs")
# plots
data = rbind(df1, df2, df3)

box_values_group(data, "Pearson Correlations Summary (Beta-ISGs Microbiome LPS)")
density_values_group(data, "Pearson Correlations Summary (Beta-ISGs Microbiome LPS)")

######## isgs ############
df1 = get_corr(mibi[-n_na, ], isgs_rlog[-n_na, ], "Microbiome:Core-ISGs")
# stats::cor(clin$LPS, mibi, use = "pairwise.complete.obs", method = "pearson")
df2 = get_corr(LPS[-n_na], mibi[-n_na, ], "LPS:Microbiome")
df3 = get_corr(LPS[-n_na], isgs_rlog[-n_na, ], "LPS:Core-ISGs")
# plots
data = rbind(df1, df2, df3)

box_values_group(data, "Pearson Correlations Summary (Core-ISGs Microbiome LPS)")
density_values_group(data, "Pearson Correlations Summmary (Core-ISGs Microbiome LPS)")


