######################
## Set up workspace
######################
rm(list = ls())
library(openxlsx)
library(tidyverse)
library(magrittr)

library(caret)   # for cross validation

library(SmCCNet)
library(parallel)

# this is K clusters
Kc = 6
cl <- makeCluster(Kc, type = "FORK") # Create K parallel threads.

options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()

# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"

k_fold_

CVDir <- "Example3foldCV/"
# create dir under getwd() 
dir.create(CVDir)
