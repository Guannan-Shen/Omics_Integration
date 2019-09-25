options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()
library(tidyverse)
library(magrittr)

# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
source( paste0(dir, "Code/CV_caret_smcca.R") )
source( paste0(dir, "Code/put_together.R") )

###############################################################
###### combination of prev_ RA cutoff 30% 40% 1% 2%
###############################################################

######33 for global transcriptome 100, 50 cut outlier fdr 0.05, 0.7, 0.9, 1000, micro genus 40% 2% ####
## sCD14 0.20 0.75
## HIV status 0.50 0.85
## LTA 0.05 0.25 (L1 is for transcriptome L2 is for microbiome)
## LPS 0.10 0.15

## CRP 0.10 0.05
## TWO OMICS 0.50 0.10

print("Transcriptome 1894, Genus Microbiome 44 (40, 2), 49 (30, 2)")

##### loading data ######
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

#### phenotype contains ID ####
clin <- rescaled_cli() 
CD14 <- clin %>% dplyr::select(CD14)
# disease status
HIV <- ifelse(clin$Group == "hiv", 1, 0) %>% as.data.frame()
colnames(HIV) <- "HIV"
#
LTA <- clin %>% dplyr::select(LTA)
#
LPS <- clin %>% select(LPS)
n_na <- which(is.na(LPS))

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
prev <- 30
ra <- 2
omcis_name <- "30_2_Global_100_50_Genus"
  
micro_data <- load_filtered_micro_level_samples("genus",  
                                                prevalence = prev, RA = ra, wd = "Ubuntu")
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(sum(mibi_outlier))

###################### get l1l2 ################
## sCD14 0.20 0.75
## HIV status 0.50 0.85
## LTA 0.05 0.25 (L1 is for transcriptome L2 is for microbiome)
## LPS 0.10 0.15

check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(LTA), K = 4)
############## sCD 14  ##############

CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          ## sCD14 0.20 0.75
          pen1 = 0.25, pen2 = 0.8,
          NoTrait = FALSE,
          bytrait = FALSE)

############## LTA ##################
CVDir <- get_CVDir(Y = LTA, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(LTA), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          # ## LTA 0.05 0.25
          pen1 = 0.2, pen2 = 0.3,
          NoTrait = FALSE,
          bytrait = FALSE)

############## LPS ##################
CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(X1 = filtered_rlog[-n_na, filtered_outlier], X2 = mibi[-n_na, mibi_outlier], 
          Y = data.frame(LPS[-n_na, ]), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          # ## LPS 0.10 0.15
          pen1 = 0.2, pen2 = 0.25,
          NoTrait = FALSE,
          bytrait = FALSE)

#############3 HIV ############
CVDir <- get_CVDir(Y = HIV, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(HIV), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          ## HIV status 0.50 0.85
          pen1 = 0.6, pen2 = 0.9,
          NoTrait = FALSE,
          bytrait = TRUE)


################################ 30 1 #######################################
#############33##3 different mibi ##############
prev <- 30
ra <- 1
omcis_name <- "30_1_Global_100_50_Genus"

micro_data <- load_filtered_micro_level_samples("genus",  
                                                prevalence = prev, RA = ra, wd = "Ubuntu")
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(sum(mibi_outlier))

###################### get l1l2 ################
## sCD14 0.20 0.75
## HIV status 0.50 0.85
## LTA 0.05 0.25 (L1 is for transcriptome L2 is for microbiome)
## LPS 0.10 0.15

check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(LTA), K = 4)
############## sCD 14  ##############

CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          ## sCD14 0.20 0.75
          pen1 = 0.25, pen2 = 0.8,
          NoTrait = FALSE,
          bytrait = FALSE)

############## LTA ##################
CVDir <- get_CVDir(Y = LTA, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(LTA), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          # ## LTA 0.05 0.25
          pen1 = 0.2, pen2 = 0.3,
          NoTrait = FALSE,
          bytrait = FALSE)

############## LPS ##################
CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(X1 = filtered_rlog[-n_na, filtered_outlier], X2 = mibi[-n_na, mibi_outlier], 
          Y = data.frame(LPS[-n_na, ]), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          # ## LPS 0.10 0.15
          pen1 = 0.2, pen2 = 0.25,
          NoTrait = FALSE,
          bytrait = FALSE)

#############3 HIV ############
CVDir <- get_CVDir(Y = HIV, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(HIV), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          ## HIV status 0.50 0.85
          pen1 = 0.6, pen2 = 0.9,
          NoTrait = FALSE,
          bytrait = TRUE)



################################ 40 1 #######################################
#############33##3 different mibi ##############
prev <- 40
ra <- 1
omcis_name <- "40_1_Global_100_50_Genus"

micro_data <- load_filtered_micro_level_samples("genus",  
                                                prevalence = prev, RA = ra, wd = "Ubuntu")
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(sum(mibi_outlier))

###################### get l1l2 ################
## sCD14 0.20 0.75
## HIV status 0.50 0.85
## LTA 0.05 0.25 (L1 is for transcriptome L2 is for microbiome)
## LPS 0.10 0.15

check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(LTA), K = 4)
############## sCD 14  ##############

CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          ## sCD14 0.20 0.75
          pen1 = 0.25, pen2 = 0.8,
          NoTrait = FALSE,
          bytrait = FALSE)

############## LTA ##################
CVDir <- get_CVDir(Y = LTA, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(LTA), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          # ## LTA 0.05 0.25
          pen1 = 0.2, pen2 = 0.3,
          NoTrait = FALSE,
          bytrait = FALSE)

############## LPS ##################
CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(X1 = filtered_rlog[-n_na, filtered_outlier], X2 = mibi[-n_na, mibi_outlier], 
          Y = data.frame(LPS[-n_na, ]), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          # ## LPS 0.10 0.15
          pen1 = 0.2, pen2 = 0.25,
          NoTrait = FALSE,
          bytrait = FALSE)

#############3 HIV ############
CVDir <- get_CVDir(Y = HIV, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(HIV), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          ## HIV status 0.50 0.85
          pen1 = 0.6, pen2 = 0.9,
          NoTrait = FALSE,
          bytrait = TRUE)





