dir = "~/Documents/gitlab/Omics_Integration/"
source( paste0(dir, "Code/ref_plots.R") )
library(outliers)

# clinical data 
source( paste0(dir, "Code/6_5_clean_clinical.R") )
clin <-  rescaled_cli() 
####### need the testing datasets 

####### plots to imprically show outliers ###
# dot_groupby(clin, clin$LPS, clin$Group, 0.8, "HIV Status", "LPS")
# dot_groupby(clin, clin$LPS, 1, 0.8, "The whole cohort", "LPS")
# 
# dot_groupby(clin, clin$CD14, clin$Group, 0.8, "HIV Status", "CD14")
# dot_groupby(clin, clin$CD14, 1, 0.8, "The whole cohort", "CD14")

## grubbs wrapper ##
grubbs_wrapper <- function(x, type, hiv){
  # print("Requires the clin dataset, group = 'hiv' or 'control' ")
  # print("Choose type from 10, 11, 20")
  library(outliers)
  if(missing(hiv)){
    grubbs = grubbs.test(x = x, 
                         type = type)
    return( ifelse(grubbs$p.value < 0.05,
                   "Outlier",
                   "Nonoutlier") )
  }
  else{
    grubbs = grubbs.test(x = x[clin$Group == hiv], 
                         type = type)
    return( ifelse(grubbs$p.value < 0.05,
                   "Outlier",
                   "Nonoutlier") )
  }
}
# ## test of functions 
# apply(isgs_rlog, 2, FUN = function(x){
#   grubbs_wrapper(x, 10)
# })
# apply(isgs_rlog, 2, var)
# apply(mibi, 2, FUN = function(x){
#   grubbs_wrapper(x, 10)
# })
# 
# sapply(colnames(mibi), FUN = function(x){
#   dot_groupby(mibi, mibi[,x], clin$Group, 0.8, "HIV Status", 
#               paste("Rescaled-clr of", x) )
# })
# apply(mibi, 2, FUN = function(x){
#   grubbs_wrapper(x, 10, "control")
# })
# ########## test of outliers ##########
# library(outliers)
# grubbs.test(clin$LPS, # a numeric vector for data values.
#             type = 10, # Integer value indicating test variant. 
#             # 10 is a test for one outlier (side is detected automatically and can be reversed by opposite parameter). 
#             # 11 is a test for two outliers on opposite tails, 20 is test for two outliers in one tail.
#             opposite = FALSE, # a logical indicating whether you want to check not the value with largest difference from the mean
#             two.sided = FALSE)
# # one outlier at each tail
# grubbs.test(clin$CD14, 
#             type = 11, 
#             opposite = FALSE, 
#             two.sided = FALSE)
# # one outlier at a side 
# grubbs.test(clin$CD14, 
#             type = 10)
# # two outliers at one side 
# grubbs.test(clin$CD14, 
#             type = 20, 
#             opposite = FALSE, 
#             two.sided = FALSE)
# ######## test outlier by group 
# grubbs.test(clin[clin$Group == 'hiv', colnames(clin) == "LPS"], 
#             type = 10)
# 
# grubbs.test(clin[clin$Group == 'control', colnames(clin) == "LPS"], 
#             type = 10)
# 
# # Test the suspicious one 
# grubbs.test(mibi$Bacte.Bacte.Prevotella[clin$Group == 'control'], 
#             type = 10)
# 
# dot_groupby(mibi, mibi$Bacte.Bacte.Prevotella, 1, 0.8, "The whole cohort", "Prevotella")
# dot_groupby(mibi, mibi$Bacte.Bacte.Prevotella, clin$Group, 0.8, "HIV Status", "Prevotella")
# length(mibi$Bacte.Bacte.Prevotella)
