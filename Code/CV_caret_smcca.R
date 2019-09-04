######################
## Set up workspace
######################
library(openxlsx)
library(tidyverse)
library(magrittr)

library(caret)   # for cross validation

library(SmCCNet)
library(parallel)

options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()

# # directory, Ubuntu 
# dir = "~/Documents/gitlab/Omics_Integration/"

############# Three testing datasets from 6_11_CV_LPS_targeted #############
###### 8_5_testing_dataset ############
## isgs_rlog, mibi, LPS

# In the code below, we show how to create CV
# data sets and check if all data sets are valid (i.e. standardizable).
check_standardize <- function(X1, X2, Y, K){
  # sample size
  n = nrow(Y)
  set.seed(12345) # Set random seed.
  foldIdx <- split(1:n, sample(1:n, K))
  for(i in 1:K){
    iIdx <- foldIdx[[i]]
    x1.train <- scale(X1[-iIdx, ])
    x2.train <- scale(X2[-iIdx, ])
    yy.train <- scale(Y[-iIdx, ])
    x1.test <- scale(X1[iIdx, ])
    x2.test <- scale(X2[iIdx, ])
    yy.test <- scale(Y[iIdx, ])
    # Check if standardized data sets are valid.
    if(is.na(min(min(x1.train), min(x2.train), min(yy.train), min(x1.test),
                 min(x2.test), min(yy.test))))
      {
      stop("Invalid scaled data. At least one of the data matrices include a
column with zero variance.")
     } else {print("We are good to go!")}
  }
}

# k_fold_lambda(X1 = isgs_rlog, X2 = mibi, Y = LPS, K = 4, CCcoef = NULL, Omics_name = "Core_ISGs_Genus")

##################################

######### get the dir name ###########
get_CVDir <- function(X1, X2, Y, K, CCcoef, Omics_name, ntrys){
  # use weights in the name
  print("K fold cross-validation!")
  clin_name = as.character(substitute(Y)) 
  name = paste0(clin_name, "_", paste(CCcoef, collapse = "") )
  # Set a CV directory.
  setwd("~/Documents/gitlab/Omics_Integration/DataProcessed/")
  CVDir <-  paste0(name, Omics_name, "_", ntrys,"_", K, "foldCV/")
  dir.create(CVDir)
  # 
  return(CVDir)
}

######### run CV ##################
CV_lambda <- function(X1, X2, Y, K, CCcoef, s1, s2, pen1, pen2){
  #
  print("K fold and K threads!")
  # Set a CV directory.
  setwd("~/Documents/gitlab/Omics_Integration/DataProcessed/")
# parameters unchanged
p1 = ncol(X1)
p2 = ncol(X2)
n = nrow(Y)
AbarLabel = c(colnames(cbind(X1, X2)))
#
print("For sample size < 30, K fold = 4 is recommended.")
print("CCcoef can be NULL or length 3 numeric vectors.")
# 
# s1 = 0.6; s2 = 0.9 # Feature sampling proportions.
SubsamplingNum = 1000 # Number of subsamples.
# Create sparsity penalty options.
pen1 = seq(.05, pen1, by = .05)
pen2 = seq(.05, pen2, by = .05)
P1P2 = expand.grid(pen1, pen2)
#
# Map (l1, l2) to (c1, c2).
c1 = sqrt(p1 * s1) * P1P2[ , 1]; c1[c1 < 1] = 1
c2 = sqrt(p2 * s2) * P1P2[ , 2]; c2[c2 < 1] = 1
# # get the name
# clin_name = as.character(substitute(Y)) 
# name = paste0(clin_name, paste(CCcoef, collapse = "") )

# CVDir <-  paste0(name, Omics_name, K, "foldCV/")
# dir.create(CVDir)
# The standardized training
# and test data sets will be saved under the CV directory.
set.seed(sample.int(1e5, 1)) # Set random seed.

foldIdx <- split(1:n, sample(1:n, K))
for(i in 1:K){
  iIdx <- foldIdx[[i]]
  x1.train <- scale(X1[-iIdx, ])
  x2.train <- scale(X2[-iIdx, ])
  yy.train <- scale(Y[-iIdx, ])
  x1.test <- scale(X1[iIdx, ])
  x2.test <- scale(X2[iIdx, ])
  yy.test <- scale(Y[iIdx, ])
  # Check if standardized data sets are valid.
  if(is.na(min(min(x1.train), min(x2.train), min(yy.train), min(x1.test),
               min(x2.test), min(yy.test)))){
    stop("Invalid scaled data. At least one of the data matrices include a
         column with zero variance.")
  }
  # Create Subdirectory
  subD <- paste0(CVDir, "CV_", i, "/")
  dir.create(subD)
  save(x1.train, x2.train, yy.train, x1.test, x2.test, yy.test,
       s1, s2, P1P2, p1, p2, SubsamplingNum, CCcoef,
       file = paste0(subD, "Data.RData"))
  }
# parallel to run cross CV
cl <- makeCluster(K, type = "FORK") # Create K parallel threads.
clusterExport(cl = cl, "CVDir") # Pass on variable CVDir to each thread.
parSapply(cl, 1:K, function(CVidx){
  # Reload source code files for each thread
  library(SmCCNet)
  # Create a result directory for each thread.
  subD <- paste0(CVDir, "CV_", CVidx, "/")
  load(paste0(subD, "Data.RData"))
  dir.create(paste0(subD, "SmCCA/"))
  
  RhoTrain <- RhoTest <- DeltaCor <- rep(0, nrow(P1P2))
  for(idx in 1:nrow(P1P2)){
    # Consider one pair of sparsity penalties at a time.
    l1 <- P1P2[idx, 1]
    l2 <- P1P2[idx, 2]
    
    # Run SmCCA on the subsamples (Figure 1, Step II)
    Ws <- getRobustPseudoWeights(x1.train, x2.train, yy.train, l1, l2, 
                                 s1, s2, NoTrait = FALSE,
                                 FilterByTrait = FALSE, 
                                 SubsamplingNum = SubsamplingNum, 
                                 CCcoef = CCcoef)
    
    # Aggregate pseudo-canonical weights from the subsamples.
    meanW <- rowMeans(Ws)
    v <- meanW[1:p1]
    u <- meanW[p1 + 1:p2]
    
    # Compute the prediction error for given CV fold and sparsity penalties.
    if(is.null(CCcoef)){CCcoef <- rep(1, 3)} # Unweighted SmCCA.
    rho.train <- cor(x1.train %*% v, x2.train %*% u) * CCcoef[1] + 
      cor(x1.train %*% v, yy.train) * CCcoef[2] + 
      cor(x2.train %*% u, yy.train) * CCcoef[3]
    rho.test <- cor(x1.test %*% v, x2.test %*% u) * CCcoef[1] +
      cor(x1.test %*% v, yy.test) * CCcoef[2] + 
      cor(x2.test %*% u, yy.test) * CCcoef[3]
    RhoTrain[idx] <- round(rho.train, digits = 5)
    RhoTest[idx] <- round(rho.test, digits = 5)
    DeltaCor[idx] <- abs(rho.train - rho.test)
    
    # Periodically save results in a temporary file.
    if(idx %% 10 == 0){
      save(P1P2, RhoTrain, RhoTest, DeltaCor, idx, 
           file = paste0(subD, "temp.RData"))
    }
  }
  # Record prediction errors for given CV fold and all sparsity penalty 
  # options.
  DeltaCor.all <- cbind(P1P2, RhoTrain, RhoTest, DeltaCor)
  colnames(DeltaCor.all) <- c("l1", "l2", "Training CC", "Test CC", 
                              "CC Pred. Error")
  write.csv(DeltaCor.all, 
            file = paste0(subD, "SmCCA/PredictionError.csv"))
  
  # Remove the temporary file.
  system(paste0("rm ", subD, "temp.RData"))
  return(CVidx)
})

# Close cluster
stopCluster(cl)
#############################################
# Combine prediction errors from all K folds and compute the total prediction
# error for each sparsity penalty pair.
testCC <- predError <- NULL
for(j in 1:K){
  resultT <- paste0(CVDir, "CV_", j, "/SmCCA/PredictionError.csv")
  dCorT <- read.csv(resultT)[ , -1]
  testCC <- cbind(testCC, abs(dCorT[ , 4]))
  predError <- cbind(predError, dCorT[ , 5])
}
S1 <- rowMeans(testCC)
S2 <- rowMeans(predError)
T12 <- dCorT[ , -3]; T12[ , 3] <- S1; T12[ , 4] <- S2
write.csv(T12, file = paste0(CVDir, "TotalPredictionError.csv"))


library(plotly)
library(reshape2)

f1 <- list(
  family = "Arial, sans-serif",
  size = 20,
  color = "black"
)
f2 <- list(
  family = "Old Standard TT, serif",
  size = 20,
  color = "black"
)
a <- list(
  title = paste("L1", "(s1 =", s1, "s2 =", s2,")",sep = " "),
  titlefont = f1,
  showticklabels = TRUE,
  tickfont = f2
)
b <- list(
  title = "L2",
  titlefont = f1,
  showticklabels = TRUE,
  tickfont = f2
)
hmelt <- melt(T12[ , -3], id.vars = c("l1", "l2"))
contourPlot <- plot_ly(hmelt, x = ~l1, y = ~l2, z = ~value, type = "contour") %>%
  layout(xaxis = a, yaxis = b, showlegend = TRUE, legend = f1)  
# orca not works for me
 export(contourPlot, paste0(CVDir, "L1L2.png"))

}

# orca(contourPlot, file = paste0(CVDir, "TotalPredictionError.pdf"))

# ####### test run ################
# check_standardize(isgs_rlog[-n_na, ], mibi[-n_na, ], LPS, K = 4)
# check_standardize(genesbeta_rlog[-n_na, ], mibi[-n_na, ], LPS, K = 4)
# 
# # CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, Omics_name = "Core_ISGs_Genus", ntrys = 1)
# # CV_lambda(X1 = isgs_rlog[-n_na, ], X2 = mibi[-n_na, ], Y = LPS, K = 4, 
# #           CCcoef = NULL, s1 = 0.9, s2= 0.9,
# #           pen1 = 0.6, pen2 = 0.3)
# 
# ################ grid search ################
# #### unweighted isgs
# s1 = seq(.5, 0.9, by = .1)
# s2 = seq(.5, 0.9, by = .1)
# s1s2 = rbind(expand.grid(s1, s2), expand.grid(s1, s2), expand.grid(s1, s2))
# 
# for (i in 1:nrow(s1s2)){
#   CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, Omics_name = "Core_ISGs_Genus", ntrys = i)
#   CV_lambda(X1 = isgs_rlog[-n_na, ], X2 = mibi[-n_na, ], Y = LPS, K = 4, 
#             CCcoef = NULL, s1 = s1s2[i,1], s2= s1s2[i,2],
#             pen1 = 0.6, pen2 = 0.3)
# }
# ### weighted isgs 
# s1 = seq(.6, 0.9, by = .1)
# s2 = seq(.7, 0.9, by = .1)
# s1s2 = rbind(expand.grid(s1, s2), expand.grid(s1, s2), expand.grid(s1, s2))
# for (i in 1:nrow(s1s2)){
#   CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = c(20,1,20), Omics_name = "Core_ISGs_Genus", ntrys = i)
#   CV_lambda(X1 = isgs_rlog[-n_na, ], X2 = mibi[-n_na, ], Y = LPS, K = 4, 
#             CCcoef = c(20,1,20), s1 = s1s2[i,1], s2= s1s2[i,2],
#             pen1 = 0.5, pen2 = 0.3)
# }
# 
# ######## updated grid search of genesbeta (narrow down) ######
# ### unweighted genesbeta
# s1 = seq(.6, 0.9, by = .1)
# s2 = seq(.7, 0.9, by = .1)
# s1s2 = rbind(expand.grid(s1, s2), expand.grid(s1, s2), expand.grid(s1, s2))
# for (i in 1:nrow(s1s2)){
#   CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = NULL, Omics_name = "genesbeta_Genus", ntrys = i)
#   CV_lambda(X1 = genesbeta_rlog[-n_na, ], X2 = mibi[-n_na, ], Y = LPS, K = 4, 
#             CCcoef = NULL, s1 = s1s2[i,1], s2= s1s2[i,2],
#             pen1 = 0.6, pen2 = 0.3)
# }
# 
# ### weighted genesbeta
# s1 = seq(.6, 0.9, by = .1)
# s2 = seq(.7, 0.9, by = .1)
# s1s2 = rbind(expand.grid(s1, s2), expand.grid(s1, s2), expand.grid(s1, s2))
# for (i in 1:nrow(s1s2)){
#   CVDir <- get_CVDir(Y = LPS, K = 4, CCcoef = c(10,1,10), Omics_name = "genesbeta_Genus", ntrys = i)
#   CV_lambda(X1 = genesbeta_rlog[-n_na, ], X2 = mibi[-n_na, ], Y = LPS, K = 4, 
#             CCcoef = c(10,1,10), s1 = s1s2[i,1], s2= s1s2[i,2],
#             pen1 = 0.6, pen2 = 0.3)
# }
