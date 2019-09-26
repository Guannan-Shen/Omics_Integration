options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()
library(tidyverse)
library(magrittr)

# directory, Ubuntu this is directory to source files 
dir = "~/Documents/gitlab/Omics_Integration/"
source( paste0(dir, "Code/put_together.R") )

############# load in clinical data ###############
# load all datasets 
# source( paste0(dir, "Code/8_5_testing_dataset.R") )

###### focus on HIV status, sCD14, LPS, LTA, crp ########3 
###### information for table 1 ##########
# at the 5_21_2019 .Rmd file
print("Table 1 Done! (5_21_2019 .Rmd)")

########## t-test and linear regression against HIV ################
print("9_24_done!")

################### features in modules (non-reduced) (correlation, association) against the phenotype #####
all_modules <- function(modules){
  n_networks = length(modules)
  sum_nodes = NULL
  # get all feature index in the 
  for (i in 1:n_networks){
    sum_nodes = c(sum_nodes, modules[[i]])
  }
  return(sum_nodes)
}

module_0 <- function(abar, modules){
  print("All features excluded in identified modules are grouped as Module 0.")
  all_nodes = all_modules(modules)
  # all features indices
  features = 1:ncol(abar)
  return(features[-all_nodes])
}

######## get new set of modules including 0 modules #########
add_module_0 <- function(abar, modules){
  n_networks = length(modules)
  mo_0 = module_0(abar, modules)
  res = vector("list", n_networks + 1)
  res[[1]] = mo_0
  for (i in 2: (n_networks + 1) ){
    res[[i]] = modules[[i-1]]
  }
  return(res)
}


signed_sim_matrix_cut <- function(X1, X2, abar, modules, edgecut){
  print("Might include module 0 in modules_0\n
        modules is used here")
  n_networks = length(modules)
  # sign from correlation 
  x = cbind(X1, X2)
  corr = cor(x, method = "pearson")
  # store length of cut modules and cut modules 
  res = vector("list", n_networks)
  for (i in 1:n_networks) {
    # get similarity matrix 
    M = as.matrix(abar[modules[[i]], modules[[i]]])
    #signed matrix 
    M = M * sign(corr[modules[[i]], modules[[i]]])
    M.node = colnames(M)
    ### edge cut 
    M[which(abs(M) < edgecut)] = 0
    newM.node = M.node[which(apply(abs(M), 1, max) > 0)]
    res[[i]] = M[newM.node, newM.node]
  }
  return(res)
}

######### test run and comparison 
# cut_modules <- signed_sim_matrix_cut(X1 = filtered_rlog[, filtered_outlier], 
#                       X2 = mibi[, mibi_outlier],
#                       abar,
#                       modules,
#                       edgecut = 0.01)
# for (i in 1:5){
#   print(dim(cut_modules[[i]]))
# }
# 
# 
# edgecut_by(X1 = filtered_rlog[, filtered_outlier], 
#            X2 = mibi[, mibi_outlier], 0.01)

############### get levels (rescaled) values of features ##########
levels_mo_feature <- function(X1, X2, modules){
  n_networks = length(modules)
  x = cbind(X1, X2)
  d_list = vector("list", n_networks)
  for (i in 1:n_networks) {
    d_list[[i]] = x[, modules[[i]]]
  }
  return(d_list)
}

######### get levels (rescaled) values of features from colnames of 
levels_mo_feature_mat <- function(X1, X2, mats){
  n_networks = length(modules)
  x = cbind(X1, X2)
  d_list = vector("list", n_networks)
  for (i in 1:n_networks) {
    mat_index = colnames(reduced_sim[[i]])
    d_list[[i]] = x[,  mat_index]
  }
  return(d_list)
}



########### get pearson correlations against the phenotype ##########
corr_list_module <- function(fea_list, pheno, modules_index, modulelabel){
  print("No missing values")
  n_networks = length(fea_list)
  res = vector("list", n_networks)
  for (i in 1:n_networks){
    # pearson correlations
    res[[i]] = stats::cor(fea_list[[i]], pheno, method = "pearson") %>% as.data.frame() %>% stack() %>% 
      as.data.frame()  %>% dplyr::mutate(ind = paste0(modulelabel, modules_index[i]) ) %>%
      dplyr::select(values, ind)
  }
  # combine data.frames from a list into one dataframe
  df = do.call("rbind", res)
  return(df)
}

########### grouped dot plot #############
corr_modules <- function(data, groupby, y, xlab, ylab){
  p = ggplot(data, aes(x=groupby, y=y)) + 
    # dot plot
    # geom_dotplot(binaxis='y', stackdir='center', dotsize = dotsize)+
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.09) +
    scale_colour_grey(start = 0.2, end = 0.8) +
    scale_y_continuous(name = waiver(), limits = c(-1.02, 1.02), breaks = c(-1, -0.5, 0, 0.5, 1), 
                       labels = c(-1, -0.5, 0, 0.5, 1)) +
    theme_bw() +
    # mean and 2 std
    # stat_summary(fun.data="mean_sdl", fun.args = list(mult=2), 
    #              geom="crossbar", width=0.3) +
    labs(x = xlab, y = ylab) +
    coord_flip() +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(p)
}


################ test run #################
# fea_list <- levels_mo_feature(X1 = filtered_rlog[, filtered_outlier], 
#                   X2 = mibi[, mibi_outlier], 
#                   module_s_0)
# 
# fea_trimmed <- levels_mo_feature_mat(X1 = filtered_rlog[, filtered_outlier], 
#                       X2 = mibi[, mibi_outlier], 
#                       reduced_sim)
# 
# test <- corr_list_module(fea_list, CD14, 0:5, modulelabel = "module")
# 
# test1 <- corr_list_module(fea_trimmed, CD14, 1:5, modulelabel = "trimmed module")
# 
# mo_cor_test <- rbind(test, test1)
# 
# corr_modules(mo_cor_test, mo_cor_test$ind, mo_cor_test$values,  "", "Correlation")
