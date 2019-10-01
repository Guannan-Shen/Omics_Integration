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
# clinical data 
source( paste0(dir, "Code/6_5_clean_clinical.R") )
# subset global rna-seq by mean, var
source( paste0(dir, "Code/generate_genelist.R") )
# clean and transform transcriptome data, subset of genes
source( paste0(dir, "Code/6_5_clean_transcriptome.R") )
# outliers 
source( paste0(dir, "Code/outliers.R") )
# source to get the load_filtered_micro_level function to get clr of RA
source( paste0(dir, "Code/5_29_Generate_filtered_Data_Microbiome.R") )

# diagnostic plots and tables
source( paste0(dir, "Code/ref_plots.R") )

###### load in ####
#### phenotype contains ID ####
clin <- rescaled_cli() 
CD14 <- clin %>% dplyr::select(CD14)
ID <- clin %>% dplyr::select(ID)
sex <- clin %>% dplyr::select(sex)
age <- clin %>% dplyr::select(age)

# disease status
HIV <- ifelse(clin$Group == "hiv", 1, 0) %>% as.data.frame()
colnames(HIV) <- "HIV"
#
LTA <- clin %>% dplyr::select(LTA)
#
LPS <- clin %>% dplyr::select(LPS)
n_na <- which(is.na(LPS))


##### Transcriptome ######
rna_isgs <- as.data.frame(read.delim( paste0(dir, "DataRaw/hiv_infected_un/coreISG")) ) %>% 
  rescaled_rna(., rlog = T)
# the isgs data
isgs_rlog <- rna_isgs[[1]]
# names
colnames(isgs_rlog) <- rna_isgs[[2]]$Symbol
######## genes beta
rna_genesbeta <- as.data.frame(read.delim( paste0(dir, 
                                                  "DataRaw/hiv_infected_un/genesbeta")) ) %>% 
  rescaled_rna(., rlog = T)
genesbeta_rlog <- rna_genesbeta[[1]]
colnames(genesbeta_rlog) <- rna_genesbeta[[2]]$Symbol

isgs_outlier <- grubbs_df(isgs_rlog, 2, 10)$fdr > 0.05
print(sum(isgs_outlier))
genesbeta_outlier <- grubbs_df(genesbeta_rlog, 2, 10)$fdr > 0.05
print(sum(genesbeta_outlier))

#### global filtered #######3
mean_cut <- 100
var_cut <- 50
dim(get_tmm()$df)
filtered_rna <- filter_rescale_rna(mean_cut, var_cut, T)
filtered_rlog <- filtered_rna[[1]]
colnames(filtered_rlog) <- filtered_rna[[2]]$Symbol
print(ncol(filtered_rlog))

filtered_outlier <- grubbs_df(filtered_rlog, 2, 10)$fdr > 0.05
print(sum(filtered_outlier))

#############33##3 different mibi ##############
prev <- 40
ra <- 2
# omcis_name <- "30_2_Global_100_50_Genus"

micro_data <- load_filtered_micro_level_samples("genus",  
                                                prevalence = prev, RA = ra, wd = "Ubuntu")
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(sum(mibi_outlier))


######## hints of weights genesbeta ##########

df1 = get_corr(mibi[-n_na, ], genesbeta_rlog[-n_na, ], "Microbiome:Beta-ISGs")
# stats::cor(clin$LPS, mibi, use = "pairwise.complete.obs", method = "pearson")
df2 = get_corr(LPS[-n_na,], mibi[-n_na, ], "LPS:Microbiome")
df3 = get_corr(LPS[-n_na,], genesbeta_rlog[-n_na, ], "LPS:Beta-ISGs")
# plots
data = rbind(df1, df2, df3)

# box_values_group(data, "Pearson Correlations Summary (Beta-ISGs Microbiome LPS)")
# density_values_group(data, "Pearson Correlations Summary (Beta-ISGs Microbiome LPS)")

######## isgs ############
df1 = get_corr(mibi[-n_na, ], isgs_rlog[-n_na, ], "Microbiome:Core-ISGs")
# stats::cor(clin$LPS, mibi, use = "pairwise.complete.obs", method = "pearson")
df2 = get_corr(LPS[-n_na,], mibi[-n_na, ], "LPS:Microbiome")
df3 = get_corr(LPS[-n_na,], isgs_rlog[-n_na, ], "LPS:Core-ISGs")
# plots
data = rbind(df1, df2, df3)

# box_values_group(data, "Pearson Correlations Summary (Core-ISGs Microbiome LPS)")
# density_values_group(data, "Pearson Correlations Summmary (Core-ISGs Microbiome LPS)")

print("Transcriptome cutoffs mean 100, variance 50, Microbiome cutoffs prevalence 40%, RA 2%")

