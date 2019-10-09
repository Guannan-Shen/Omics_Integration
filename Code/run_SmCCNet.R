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
source( paste0(dir, "Code/wrappers.R") )
# run smccnet
source( paste0(dir, "Code/put_together.R") )

########## with CD14 #############3
setwd("~/Documents/gitlab/Omics_Integration/DataProcessed/")
CVDir <- "CD14_Outlier3_Global_100_50_Genus_3_4foldCV/"



########## without CD14 ###########
run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
            X2 = mibi[, mibi_outlier],
            Y = NULL,
            l1 = 0.5, 
            l2 = 0.1, 
            s1 = 0.7, 
            s2 = 0.9, 
            weights = NULL,
            # n_na = n_na,
            # NoTrait itself is to control whether to use Y or not 
            NoTrait = TRUE,
            EdgeCut = 0,
            bytrait = FALSE)