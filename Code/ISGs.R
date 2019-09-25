options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()
library(tidyverse)
library(magrittr)

# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
source( paste0(dir, "Code/CV_caret_smcca.R") )
source( paste0(dir, "Code/put_together.R") )

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

######### Transcriptome ###############
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

################## 40 2 microbiome genus ###################
prev <- 40
ra <- 2
omcis_name <- "40_2_core_ISGs_Genus_rerun"

micro_data <- load_filtered_micro_level_samples("genus",  
                                                prevalence = prev, RA = ra, wd = "Ubuntu")
micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(sum(mibi_outlier))


##################### isgs ###############33
###################### get l1l2 ################

check_standardize(isgs_rlog[, isgs_outlier], mibi[, mibi_outlier], 
                  data.frame(LTA), K = 4)
############## sCD 14  ##############

CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(isgs_rlog[, isgs_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.8, s2= 0.9,
          ## sCD14 0.20 0.75
          pen1 = 0.7, pen2 = 0.7,
          NoTrait = FALSE,
          bytrait = FALSE)

############## LTA ##################
CVDir <- get_CVDir(Y = LTA, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(isgs_rlog[, isgs_outlier], mibi[, mibi_outlier], 
          data.frame(LTA), 
          K = 4,
          CCcoef = NULL, s1 = 0.8, s2= 0.9,
          # ## LTA 0.05 0.25
          pen1 = 0.7, pen2 = 0.7,
          NoTrait = FALSE,
          bytrait = FALSE)

############## LPS ##################
CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(X1 = isgs_rlog[-n_na, isgs_outlier], X2 = mibi[-n_na, mibi_outlier], 
          Y = data.frame(LPS[-n_na, ]), 
          K = 4,
          CCcoef = NULL, s1 = 0.8, s2= 0.9,
          # ## LPS 0.10 0.15
          pen1 = 0.7, pen2 = 0.7,
          NoTrait = FALSE,
          bytrait = FALSE)

#############3 HIV ############
CVDir <- get_CVDir(Y = HIV, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(isgs_rlog[, isgs_outlier], mibi[, mibi_outlier], 
          data.frame(HIV), 
          K = 4,
          CCcoef = NULL, s1 = 0.8, s2= 0.9,
          ## HIV status 0.50 0.85
          pen1 = 0.8, pen2 = 0.8,
          NoTrait = FALSE,
          bytrait = TRUE)


##################### genesbeta ###############33
###################### get l1l2 ################
omcis_name <- "40_2_core_beta_Genus_rerun"

check_standardize(genesbeta_rlog[, genesbeta_outlier], mibi[, mibi_outlier], 
                  data.frame(LTA), K = 4)
############## sCD 14  ##############

CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(genesbeta_rlog[, genesbeta_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.8, s2= 0.9,
          ## sCD14 0.20 0.75
          pen1 = 0.7, pen2 = 0.7,
          NoTrait = FALSE,
          bytrait = FALSE)

############## LTA ##################
CVDir <- get_CVDir(Y = LTA, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(genesbeta_rlog[, genesbeta_outlier], mibi[, mibi_outlier], 
          data.frame(LTA), 
          K = 4,
          CCcoef = NULL, s1 = 0.8, s2= 0.9,
          # ## LTA 0.05 0.25
          pen1 = 0.7, pen2 = 0.7,
          NoTrait = FALSE,
          bytrait = FALSE)

############## LPS ##################
CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(X1 = genesbeta_rlog[-n_na, genesbeta_outlier], X2 = mibi[-n_na, mibi_outlier], 
          Y = data.frame(LPS[-n_na, ]), 
          K = 4,
          CCcoef = NULL, s1 = 0.8, s2= 0.9,
          # ## LPS 0.10 0.15
          pen1 = 0.7, pen2 = 0.7,
          NoTrait = FALSE,
          bytrait = FALSE)

#############3 HIV ############
CVDir <- get_CVDir(Y = HIV, K = 4, CCcoef = NULL, 
                   Omics_name = omcis_name, ntrys = 1)

CV_lambda(genesbeta_rlog[, genesbeta_outlier], mibi[, mibi_outlier], 
          data.frame(HIV), 
          K = 4,
          CCcoef = NULL, s1 = 0.8, s2= 0.9,
          ## HIV status 0.50 0.85
          pen1 = 0.8, pen2 = 0.8,
          NoTrait = FALSE,
          bytrait = TRUE)
