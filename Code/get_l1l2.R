# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
source( paste0(dir, "Code/CV_caret_smcca.R") )
source( paste0(dir, "Code/put_together.R") )
#####for 9_3_outliers ##############
check_standardize(isgs_rlog[-n_na, isgs_outlier], mibi[-n_na, mibi_outlier], 
                  data.frame(LPS[-n_na, ]), K = 4)

CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, Omics_name = "Outlier1_Core_ISGs_Genus", ntrys = 1)
CV_lambda(X1 = isgs_rlog[-n_na, isgs_outlier], X2 = mibi[-n_na, mibi_outlier], 
          Y = data.frame(LPS[-n_na, ]), 
          K = 4,
          CCcoef = NULL, s1 = 0.8, s2= 0.9,
          pen1 = 0.6, pen2 = 0.6,
          NoTrait = FALSE)

########### CD14 ############ 9_12_Global_rna_integration 
# For transciptome global level 0.7, 0.9
check_standardize(isgs_rlog[, isgs_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)
check_standardize(genesbeta_rlog[, genesbeta_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)
check_standardize(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
                  data.frame(CD14), K = 4)
####### CD14, genus microbiome 40, 2, global rna mean 100, var 50, outlier fdr 0.05 #######
CVDir <- get_CVDir(Y = CD14, K = 4, CCcoef = NULL, 
                   Omics_name = "Outlier3_Global_100_50_Genus", ntrys = 3)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          pen1 = 0.3, pen2 = 0.9,
          NoTrait = FALSE)

####### CRP, genus microbiome 40, 2, global rna mean 100, var 50, outlier fdr 0.05 #######
CVDir <- get_CVDir(Y = CRP, K = 4, CCcoef = NULL, 
                   Omics_name = "Outlier1_Global_100_50_Genus", ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CRP), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          pen1 = 0.8, pen2 = 0.3,
          NoTrait = FALSE)

############ LPS genus microbiome 40, 2, global rna mean 100, var 50, outlier fdr 0.05 #######
CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, 
                   Omics_name = "Outlier1_Global_100_50_Genus", ntrys = 1)

CV_lambda(X1 = filtered_rlog[-n_na, filtered_outlier], X2 = mibi[-n_na, mibi_outlier], 
          Y = data.frame(LPS[-n_na, ]), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          pen1 = 0.5, pen2 = 0.7,
          NoTrait = FALSE)

####### without CD14, genus microbiome 40, 2, global rna mean 100, var 50, outlier fdr 0.05 #######
CVDir <- get_CVDir(Y = NULL, K = 4, CCcoef = NULL, 
                   Omics_name = "Outlier1_Global_100_50_Genus", ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(CD14), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          pen1 = 0.6, pen2 = 0.6,
          NoTrait = TRUE)

###### with HIV status ###########

CVDir <- get_CVDir(Y = HIV, K = 4, CCcoef = NULL, 
                   Omics_name = "Outlier_Global_100_50_Genus", ntrys = 1)

CV_lambda(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 
          data.frame(HIV), 
          K = 4,
          CCcoef = NULL, s1 = 0.7, s2= 0.9,
          pen1 = 0.8, pen2 = 0.8,
          NoTrait = FALSE,
          bytrait = TRUE)

# ### withCD14 global
# run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
#             X2 = mibi[, mibi_outlier],
#             Y = CD14,
#             l1 = 0.4, 
#             l2 = 0.1, 
#             s1 = 0.7, 
#             s2 = 0.9, 
#             weights = NULL,
#             # n_na = n_na,
#             # NoTrait itself is to control whether to use Y or not 
#             NoTrait = FALSE,
#             EdgeCut = 0)

# ############## test run #################
# X1 = filtered_rlog[, filtered_outlier]
# X2 = mibi[, mibi_outlier]
# Y = NULL
# l1 = 0.3
# l2 = 0.1
# s1 = 0.7
# s2 = 0.9
# weights = NULL
# # n_na = n_na,
# # NoTrait itself is to control whether to use Y or not 
# NoTrait = TRUE
# EdgeCut = 0
# 
# x = cbind(X1, X2)
# corr = cor(x, method = "pearson")
# p1 = ncol(X1)
# p2 = ncol(X2)
# n = nrow(X1)
# AbarLabel = c(colnames(cbind(X1, X2)))
# set.seed(123)
# 
# Ws = getRobustPseudoWeights(X1, X2, Y, Lambda1 = l1,
#                             Lambda2 = l2, s1 = s1, s2 = s2, NoTrait = NoTrait, FilterByTrait = FALSE,
#                             SubsamplingNum = 1000, CCcoef = weights, trace = FALSE)
# 
# abar = getAbar(Ws, P1 = p1, FeatureLabel = AbarLabel)
# 
# modules = getMultiOmicsModules(abar, P1 = p1, PlotTree = T, CutHeight = 1 - 0.1^10)
