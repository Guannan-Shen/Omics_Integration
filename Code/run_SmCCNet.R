library(tidyverse)
library(magrittr)
library(SmCCNet)

options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()

# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
# small functions, %nin%, %||%, isnothing
source( paste0(dir, "Code/small_fun.R") )
# wrappers
# source( paste0(dir, "Code/wrappers.R") )
# run smccnet
source( paste0(dir, "Code/put_together.R") )

################################################################
################################ Genus #########################3
##############################################################33


setwd("~/Documents/gitlab/Omics_Integration/DataProcessed/")
top <- "~/Documents/gitlab/Omics_Integration/DataProcessed/"

print("Genus 20% 1%")

########## with CD14 #############3
CVDir <- "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% genus
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.2, 
            l2 = 0.45, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

########### Once we have ran SmCCNet #########3

# ######3 edge Cut 0.1 
# load(paste0(dir, "SmCCNetWeights.RData"))
# edgecut_by(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 0.1)

########## without CD14 Two Omics ###########
CVDir <- "_Unclassified_Genus_Global_100_50_20_1_4_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% genus
mibi[, mibi_outlier] %>% ncol()
  
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = NULL,
            l1 = 0.05, 
            l2 = 0.45, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = TRUE,
            EdgeCut = 0,
            bytrait = FALSE)


########## HIV ###########
CVDir <- "HIV_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% genus
mibi[, mibi_outlier] %>% ncol()

run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = HIV,
            l1 = 0.05, 
            l2 = 1, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = TRUE)

########## IFNb ###########
CVDir <- "IFNb_Unclassified_Genus_Global_100_50_20_1_3_3foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% genus
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = IFNb,
            l1 = 0.05, 
            l2 = 0.25, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            n_na = n_na_IFNb,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

########## LPS ###########
CVDir <- "LPS_Unclassified_Genus_Global_100_50_20_1_3_2foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% genus
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = LPS,
            l1 = 0.75, 
            l2 = 0.45, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            n_na = n_na_LPS,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

########## with LTA #############3
CVDir <- "LTA_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% genus
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = LTA,
            l1 = 0.45,
            l2 = 0.5, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

################################################################
################################ family #########################3
##############################################################33
print("Family 20% 1%, 44 Taxa")
mibi[, mibi_outlier] %>% ncol()

########## with CD14 #############3
CVDir <- "CD14_Unclassified_Family_Global_100_50_20_1_3_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.4, 
            l2 = 0.4, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

########### Once we have ran SmCCNet #########3

# ######3 edge Cut 0.1 
# load(paste0(dir, "SmCCNetWeights.RData"))
# edgecut_by(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 0.1)

########## without CD14 Two Omics ###########
CVDir <- "_Unclassified_Family_Global_100_50_20_1_4_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = NULL,
            l1 = 0.1, 
            l2 = 0.35, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = TRUE,
            EdgeCut = 0,
            bytrait = FALSE)
# paste0(dir, "SmCCNetWeights.RData")
# load(paste0(dir, "SmCCNetWeights.RData"))
# modules[[4]]

########## HIV ###########
CVDir <- "HIV_Unclassified_Family_Global_100_50_20_1_3_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family
mibi[, mibi_outlier] %>% ncol()

run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = HIV,
            l1 = 0.1, 
            l2 = 0.9, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = TRUE)

########## IFNb ###########
CVDir <- "IFNb_Unclassified_Family_Global_100_50_20_1_3_3foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = IFNb,
            l1 = 0.1, 
            l2 = 0.2, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            n_na = n_na_IFNb,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

########## LPS ###########
CVDir <- "LPS_Unclassified_Family_Global_100_50_20_1_3_2foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = LPS,
            l1 = 0.2, 
            l2 = 0.15, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            n_na = n_na_LPS,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

########## with LTA #############3
CVDir <- "LTA_Unclassified_Family_Global_100_50_20_1_3_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = LTA,
            l1 = 0.65,
            l2 = 0.05, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)




####################### Different Mibi CD14 #####################3

############ genus #############
print("Genus 20% 1%, 70 Taxa; Genus 20% 3%, 48 Taxa")

CVDir <- "CD14_Unclassified_Genus_Global_100_50_20_3_1_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.25, 
            l2 = 0.5, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

print("Genus 40% 1%, 59 Taxa")
CVDir <- "CD14_Unclassified_Genus_Global_100_50_40_1_1_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.65, 
            l2 = 0.05, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

print("Genus 40% 1%, 59 Taxa")
CVDir <- "CD14_Unclassified_Genus_Global_100_50_40_1_0.34_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.3, 
            l2 = 0.4, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)


print("Genus 40% 3%, 42 Taxa")
CVDir <- "CD14_Unclassified_Genus_Global_100_50_40_3_1_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.6, 
            l2 = 0.1, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)


print("Genus 60% 1%, 42 Taxa")
CVDir <- "CD14_Unclassified_Genus_Global_100_50_60_1_1_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.3, 
            l2 = 0.5, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)


print("Genus 60% 3%, 33 Taxa")
CVDir <- "CD14_Unclassified_Genus_Global_100_50_60_3_1_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.15, 
            l2 = 0.55, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

################################
########## family ##########

print("Family 20% 1%, 44 Taxa; Family 20% 3%, 33 Taxa")

CVDir <- "CD14_Unclassified_Family_Global_100_50_20_3_1_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.55, 
            l2 = 0.4, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

print("Family 40% 1%, 38 Taxa")
CVDir <- "CD14_Unclassified_Family_Global_100_50_40_1_1_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.3, 
            l2 = 0.7, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

print("Family 40% 3%, 27 Taxa")
CVDir <- "CD14_Unclassified_Family_Global_100_50_40_3_1_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.6, 
            l2 = 0.7, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

print("Family 60% 1%, 29 Taxa")
CVDir <- "CD14_Unclassified_Family_Global_100_50_60_1.1_1_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.7, 
            l2 = 0.8, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

print("Family 60% 3%, 23 Taxa")
CVDir <- "CD14_Unclassified_Family_Global_100_50_60_3_1_4foldCV/"
#### the results of run_SmCCNet will be saved in this folder
dir <- paste0(top, CVDir)
########## check the mibi using, 20% 1% family 
mibi[, mibi_outlier] %>% ncol()

### the wrapper of SmCCNet, run SmCCNet with L1, L2
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = CD14,
            l1 = 0.3, 
            l2 = 0.2, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na for missing phenotype values in some subjects
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = FALSE,
            EdgeCut = 0,
            bytrait = FALSE)

