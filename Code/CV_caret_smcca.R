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

options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()

# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"

k_fold_lambda <- function(X1, X2, Y, Y_name, K, K_cluster, CCcoef){
  # parameters unchanged
  p1 = ncol(X1)
  p2 = ncol(X2)
  n = nrow(X1)
  AbarLabel = c(colnames(cbind(X1, X2)))
  #
  print("For sample size < 30, K fold = 4 is recommended.")
  print("CCcoef can be NULL or length 3 numeric vectors.")
  # 
  s1 = 0.7; s2 = 0.9 # Feature sampling proportions.
  SubsamplingNum = 1000 # Number of subsamples.
  # Create sparsity penalty options.
  pen1 = seq(.05, .3, by = .05)
  pen2 = seq(.05, .3, by = .05)
  P1P2 = expand.grid(pen1, pen2)
  #
  # Map (l1, l2) to (c1, c2).
  c1 <- sqrt(p1 * s1) * P1P2[ , 1]; c1[c1 < 1] = 1
  c2 <- sqrt(p2 * s2) * P1P2[ , 2]; c2[c2 < 1] = 1
  
}


# this is K clusters
Kc = 6
cl <- makeCluster(Kc, type = "FORK") # Create K parallel threads.

setwd("~/Documents/gitlab/Omics_Integration/DataProcessed/")
CVDir <- "Example3foldCV/"
# create dir under getwd() 
dir.create(CVDir)

split(1:n, sample(1:n, K))
split(1:25, sample(1:25, 3))
