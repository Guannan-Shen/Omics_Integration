# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
source( paste0(dir, "Code/CV_caret_smcca.R") )
#####for 9_3_outliers ##############
check_standardize(isgs_rlog[-n_na, isgs_outlier], mibi[-n_na, mibi_outlier], 
                  data.frame(LPS[-n_na, ]), K = 4)

CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, Omics_name = "Outlier1_Core_ISGs_Genus", ntrys = 1)
CV_lambda(X1 = isgs_rlog[-n_na, isgs_outlier], X2 = mibi[-n_na, mibi_outlier], 
          Y = data.frame(LPS[-n_na, ]), 
          K = 4,
          CCcoef = NULL, s1 = 0.8, s2= 0.9,
          pen1 = 0.6, pen2 = 0.6)

## 