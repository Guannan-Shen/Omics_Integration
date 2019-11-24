rm(list = ls())
library(openxlsx)
library(broom)
library(tools)
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

##### Transcriptome #####
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
##### "genus" or "family"
micro_level <- "genus"
prev <- 20
ra <- 3

# omcis_name <- "30_2_Global_100_50_Genus"
micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

##################### get L1L2 #########################
# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
source( paste0(dir, "Code/CV_caret_smcca.R") )
source( paste0(dir, "Code/put_together.R") )

# For transciptome global level 0.7, 0.9
check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)

mibi[, mibi_outlier] %>% ncol() %>% print()

##### 2222 #############
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = paste("Unclassified", tools::toTitleCase(micro_level), 
                                       "Global", mean_cut, var_cut, prev, ra, sep = "_"), 
                   ntrys = 1)

######## the L1 L2 contour plot will be saved to CVDir #########
CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          c_pen1 = c(0.05, 0.9), c_pen2 = c(0.05, 0.9),
          NoTrait = FALSE,
          bytrait = FALSE)

#################### diff cutoffs ###############3
prev <- 40
ra <- 1
#################################################

# omcis_name <- "30_2_Global_100_50_Genus"
micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

# For transciptome global level 0.7, 0.9
check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)

mibi[, mibi_outlier] %>% ncol() %>% print()

##### 2222 #############
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = paste("Unclassified", tools::toTitleCase(micro_level), 
                                      "Global", mean_cut, var_cut, prev, ra, sep = "_"), 
                   ntrys = 1)

######## the L1 L2 contour plot will be saved to CVDir #########
CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          c_pen1 = c(0.05, 0.9), c_pen2 = c(0.05, 0.9),
          NoTrait = FALSE,
          bytrait = FALSE)


#################### diff cutoffs ###############3
prev <- 40
ra <- 3
#################################################

# omcis_name <- "30_2_Global_100_50_Genus"
micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

# For transciptome global level 0.7, 0.9
check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)

mibi[, mibi_outlier] %>% ncol() %>% print()

##### 2222 #############
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = paste("Unclassified", tools::toTitleCase(micro_level), 
                                      "Global", mean_cut, var_cut, prev, ra, sep = "_"), 
                   ntrys = 1)

######## the L1 L2 contour plot will be saved to CVDir #########
CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          c_pen1 = c(0.05, 0.9), c_pen2 = c(0.05, 0.9),
          NoTrait = FALSE,
          bytrait = FALSE)

#################### diff cutoffs ###############3
prev <- 60
ra <- 1
#################################################

# omcis_name <- "30_2_Global_100_50_Genus"
micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

# For transciptome global level 0.7, 0.9
check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)

mibi[, mibi_outlier] %>% ncol() %>% print()

##### 2222 #############
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = paste("Unclassified", tools::toTitleCase(micro_level), 
                                      "Global", mean_cut, var_cut, prev, ra, sep = "_"), 
                   ntrys = 1)

######## the L1 L2 contour plot will be saved to CVDir #########
CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          c_pen1 = c(0.05, 0.9), c_pen2 = c(0.05, 0.9),
          NoTrait = FALSE,
          bytrait = FALSE)


#################### diff cutoffs ###############3
prev <- 60
ra <- 3
#################################################

# omcis_name <- "30_2_Global_100_50_Genus"
micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

# For transciptome global level 0.7, 0.9
check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)

mibi[, mibi_outlier] %>% ncol() %>% print()

##### 2222 #############
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = paste("Unclassified", tools::toTitleCase(micro_level), 
                                      "Global", mean_cut, var_cut, prev, ra, sep = "_"), 
                   ntrys = 1)

######## the L1 L2 contour plot will be saved to CVDir #########
CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          c_pen1 = c(0.05, 0.9), c_pen2 = c(0.05, 0.9),
          NoTrait = FALSE,
          bytrait = FALSE)


######################3 for family ########################
micro_level <- "family"
prev <- 20
ra <- 3

# omcis_name <- "30_2_Global_100_50_Genus"
micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

##### 2222 #############
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = paste("Unclassified", tools::toTitleCase(micro_level), 
                                      "Global", mean_cut, var_cut, prev, ra, sep = "_"), 
                   ntrys = 1)

######## the L1 L2 contour plot will be saved to CVDir #########
CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          c_pen1 = c(0.05, 0.9), c_pen2 = c(0.05, 0.9),
          NoTrait = FALSE,
          bytrait = FALSE)


#################### diff cutoffs ###############3
prev <- 40
ra <- 1
#################################################

# omcis_name <- "30_2_Global_100_50_Genus"
micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

# For transciptome global level 0.7, 0.9
check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)

mibi[, mibi_outlier] %>% ncol() %>% print()

##### 2222 #############
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = paste("Unclassified", tools::toTitleCase(micro_level), 
                                      "Global", mean_cut, var_cut, prev, ra, sep = "_"), 
                   ntrys = 1)

######## the L1 L2 contour plot will be saved to CVDir #########
CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          c_pen1 = c(0.05, 0.9), c_pen2 = c(0.05, 0.9),
          NoTrait = FALSE,
          bytrait = FALSE)


#################### diff cutoffs ###############3
prev <- 40
ra <- 3
#################################################

# omcis_name <- "30_2_Global_100_50_Genus"
micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

# For transciptome global level 0.7, 0.9
check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)

mibi[, mibi_outlier] %>% ncol() %>% print()

##### 2222 #############
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = paste("Unclassified", tools::toTitleCase(micro_level), 
                                      "Global", mean_cut, var_cut, prev, ra, sep = "_"), 
                   ntrys = 1)

######## the L1 L2 contour plot will be saved to CVDir #########
CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          c_pen1 = c(0.05, 0.9), c_pen2 = c(0.05, 0.9),
          NoTrait = FALSE,
          bytrait = FALSE)

#################### diff cutoffs ###############3
prev <- 60
ra <- 1.1
#################################################

# omcis_name <- "30_2_Global_100_50_Genus"
micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

# For transciptome global level 0.7, 0.9
check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)

mibi[, mibi_outlier] %>% ncol() %>% print()

##### 2222 #############
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = paste("Unclassified", tools::toTitleCase(micro_level), 
                                      "Global", mean_cut, var_cut, prev, ra, sep = "_"), 
                   ntrys = 1)

######## the L1 L2 contour plot will be saved to CVDir #########
CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          c_pen1 = c(0.05, 0.9), c_pen2 = c(0.05, 0.9),
          NoTrait = FALSE,
          bytrait = FALSE)


#################### diff cutoffs ###############3
prev <- 60
ra <- 3
#################################################

# omcis_name <- "30_2_Global_100_50_Genus"
micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])

# For transciptome global level 0.7, 0.9
check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)

mibi[, mibi_outlier] %>% ncol() %>% print()

##### 2222 #############
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = paste("Unclassified", tools::toTitleCase(micro_level), 
                                      "Global", mean_cut, var_cut, prev, ra, sep = "_"), 
                   ntrys = 1)

######## the L1 L2 contour plot will be saved to CVDir #########
CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          c_pen1 = c(0.05, 0.9), c_pen2 = c(0.05, 0.9),
          NoTrait = FALSE,
          bytrait = FALSE)


