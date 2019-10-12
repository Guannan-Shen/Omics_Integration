rm(list = ls())
library(openxlsx)
library(broom)
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
n_na_LPS <- which(is.na(LPS))
## IFNb
IFNb <- clin %>% dplyr::select(IFNb)
n_na_IFNb <-  which(is.na(IFNb))

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
print( colnames (filtered_rlog) [!filtered_outlier])
#############33##3 different mibi ##############
prev <- 20
ra <- 1
##### "genus" or "family"
micro_level <- "genus"
# omcis_name <- "30_2_Global_100_50_Genus"

micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu")
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

a <- colnames (mibi)
############# 20, 1 ############
prev <- 20
ra <- 0
micro_level <- "genus"
# omcis_name <- "30_2_Global_100_50_Genus"

micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu")
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

colnames (mibi) [colnames (mibi) %nin% a]
###### summary ############
print(paste0("Transcriptome cutoffs mean: ", mean_cut, " variance: ", var_cut, 
             " Microbiome cutoffs prevalence: ", prev, " RA: ", ra, " ", micro_level ) )
mibi[, "Other"]

mibi <- 
  
  mibi %>% dplyr::select(-Other) %>% dim()
################### outliers ################
# df <- micro_data[[1]][, "Prote:Epsil:Helicobacter"] %>% as.data.frame() %>% 
#                     dplyr::mutate(Taxa = " Proteobacteria:Epsilonproteobacteria:Helicobacter")  %>%
#                          set_colnames(c("RA", "Taxa"))
# dot_group_violin(df, df$Taxa, df$RA, "Taxa", "RA")
# dot_group_box(df, df$Taxa, df$RA, "Taxa", "RA")

# df <- filtered_rlog[, !filtered_outlier] %>% stack() %>% as.data.frame() 
# dot_group_box(df, df$ind, df$values, "Genes", "Expression Levels (rlog)")

# hist(mibi[,"Prote.Epsil.Helicobacter"])
# hist(micro_data[[1]][,"Prote:Epsil:Helicobacter"])
# micro_data[[1]] %>% colnames()
# hist(micro_data[[1]][,"Bacte:Bacte:Alistipes"])

######## hints of weights genesbeta ##########
# df5 = get_corr_1(filtered_rlog[, filtered_outlier], "Within-Transcriptome")
# summary(df5$values)
# df5[df5$values == 1, ]
# df4 = get_corr_1(mibi[, mibi_outlier], "Within-Microbiome")
# summary(df4$values)
# df4[df4$values == 1, ]

# compare_corr(mibi[, mibi_outlier], filtered_rlog[, filtered_outlier],
#              CD14, "Genus", "sCD14")
# compare_corr(mibi[, mibi_outlier], filtered_rlog[, filtered_outlier],
#              LTA, "Genus", "LTA")
# compare_corr(mibi[-n_na_LPS, mibi_outlier], 
#              filtered_rlog[-n_na_LPS, filtered_outlier],
#              LPS[-n_na_LPS, ], "Genus", "LPS")
# compare_corr(mibi[-n_na_IFNb, mibi_outlier], 
#              filtered_rlog[-n_na_IFNb, filtered_outlier],
#              IFNb[-n_na_IFNb, ], "Genus", "IFN-Beta")

# # stats::cor(clin$LPS, mibi, use = "pairwise.complete.obs", method = "pearson")

# 
# ######## isgs ############
# df1 = get_corr(mibi[-n_na, ], isgs_rlog[-n_na, ], "Microbiome:Core-ISGs")
# # stats::cor(clin$LPS, mibi, use = "pairwise.complete.obs", method = "pearson")
# df2 = get_corr(LPS[-n_na,], mibi[-n_na, ], "LPS:Microbiome")
# df3 = get_corr(LPS[-n_na,], isgs_rlog[-n_na, ], "LPS:Core-ISGs")
# # plots
# data = rbind(df1, df2, df3)

# box_values_group(data, "Pearson Correlations Summary (Core-ISGs Microbiome LPS)")
# density_values_group(data, "Pearson Correlations Summmary (Core-ISGs Microbiome LPS)")




########### select mutate ends_with ###########
###### combine unclassified to Other , at the Relative abundance level
######## check rowSums to 1, by library (subject)
##### then get the clr 

# df <- data.frame(matrix(1, ncol = 7, nrow = 4))
# colnames(df) <- c("a", "aceae", "ales", "Cyano.4C0d.2", "Bacte.Bacte.S24.7", 
#                   "Prote.Betaproteobacteria", "other")
# remove_non <- colnames(df)[4:6]
# df %>% dplyr::mutate(other = other + rowSums(df[, colnames(df)[4:6]]) + 
#                        rowSums(dplyr::select(df, ends_with("aceae")) ) +
#                        rowSums(dplyr::select(df, ends_with("ales")) ) ) %>% 
#                   dplyr::select(-remove_non, -ends_with("aceae"), -ends_with("ales"))


# G_40_2_54 <- colnames(mibi)
# 
# # using ! to not, reverse true false
# G_40_2_out <- colnames(mibi[, !mibi_outlier])
# 
# ############# or other cutoffs #########333
# prev <- 30
# ra <- 1
# # omcis_name <- "30_2_Global_100_50_Genus"
# 
# micro_data <- load_filtered_micro_level_samples("genus",  
#                                                 prevalence = prev, RA = ra, wd = "Ubuntu")
# micro_clr <- micro_data[[2]] %>% as.data.frame()
# # rescale to mean 0 and variance 1 
# mibi <- rescale_microbiome(micro_clr)
# print(ncol(mibi))
# mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
# print(sum(mibi_outlier))
# colnames(mibi)[colnames(mibi) %nin% G_40_2_54 ]
# G_30_1_66 <- colnames(mibi)
# 
# prev <- 20
# ra <- 1
# # omcis_name <- "30_2_Global_100_50_Genus"
# 
# micro_data <- load_filtered_micro_level_samples("genus",  
#                                                 prevalence = prev, RA = ra, wd = "Ubuntu")
# micro_clr <- micro_data[[2]] %>% as.data.frame()
# # rescale to mean 0 and variance 1 
# mibi <- rescale_microbiome(micro_clr)
# print(ncol(mibi))
# mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
# print(sum(mibi_outlier))
# colnames(mibi)[colnames(mibi) %nin% G_30_1_66 ]
