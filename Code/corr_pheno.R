options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()
library(tidyverse)
library(magrittr)

# directory, Ubuntu this is directory to source files 
dir = "~/Documents/gitlab/Omics_Integration/"
source( paste0(dir, "Code/put_together.R") )
source( paste0(dir, "Code/pca_getPC1.R") )

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
  print("All features not in identified modules are grouped as Module 0.")
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
  n_networks = length(mats)
  x = cbind(X1, X2)
  d_list = vector("list", n_networks)
  for (i in 1:n_networks) {
    mat_index = colnames(mats[[i]])
    d_list[[i]] = x[,  mat_index]
  }
  return(d_list)
}

# stats::cor(fea_list[[1]], LPS, method = "pearson"
#            , use = "pairwise.complete.obs" 
# )

########### get pearson correlations against the phenotype ##########
corr_list_module <- function(fea_list, pheno, modules_index, modulelabel){
  print("No missing values")
  n_networks = length(fea_list)
  res = vector("list", n_networks)
  for (i in 1:n_networks){
    # pearson correlations
    res[[i]] = stats::cor(fea_list[[i]], pheno, method = "pearson"
                          
                          ) %>% as.data.frame() %>% stack() %>% 
      as.data.frame()  %>% dplyr::mutate(ind = paste0(modulelabel, modules_index[i]) ) %>%
      dplyr::select(values, ind)
  }
  # combine data.frames from a list into one dataframe
  df = do.call("rbind", res)
  return(df)
}

corr_list_module_noNA <- function(fea_list, pheno, modules_index, modulelabel, n_na){
  print("No missing values")
  n_networks = length(fea_list)
  res = vector("list", n_networks)
  for (i in 1:n_networks){
    # pearson correlations
    res[[i]] = stats::cor(fea_list[[i]][-n_na,], pheno[-n_na,], method = "pearson"
                          
    ) %>% as.data.frame() %>% stack() %>% 
      as.data.frame()  %>% dplyr::mutate(ind = paste0(modulelabel, modules_index[i]) ) %>%
      dplyr::select(values, ind)
  }
  # combine data.frames from a list into one dataframe
  df = do.call("rbind", res)
  return(df)
}

corr_list_module_NA <- function(fea_list, pheno, modules_index, modulelabel){
  print("No missing values")
  n_networks = length(fea_list)
  res = vector("list", n_networks)
  for (i in 1:n_networks){
    # pearson correlations
    res[[i]] = stats::cor(fea_list[[i]], pheno, method = "pearson"
                          , use = "pairwise.complete.obs" 
    ) %>% as.data.frame() %>% stack() %>% 
      as.data.frame()  %>% dplyr::mutate(ind = paste0(modulelabel, modules_index[i]) ) %>%
      dplyr::select(values, ind)
  }
  # combine data.frames from a list into one dataframe
  df = do.call("rbind", res)
  return(df)
}

########### grouped dot plot #############
corr_modules <- function(data, groupby, y, xlab, ylab, run, dir){
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
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18)) +
    theme(text=element_text(family="Arial"))
  print(p)
  ggsave(filename = paste0(run, "corre_with_pheno.tiff"),device = NULL,
         path = dir, dpi = 300, compression = "lzw" )
}

########### partial correlation ##############
residual_df <- function(df, adjustfor){
    mat = apply(df, 2, function(x){
        residuals(lm(x ~ adjustfor ))
      })
    n_obs = nrow(df)
    ######## merge list of diff length to dataframe #######3
    seq.max = seq_len(max(n_obs))
    res =  sapply(mat, "[", i = seq.max)
    return(res)
  }

############## summary #############
corr_against_Y <- function(abar, modules, X1, X2, Y, run, n_strong_modules, dir, edges, Y_name){
  print("Please load the abar, modules and Ws first, this one needs the Y!")
  print(run)
  print("Check the dimension of microbiome dataset: ")
  X2 %>% ncol() %>% print()
  # number of feature
  p1 <-  ncol(X1)
  p2 <-  ncol(X2)
  n_networks <- length(modules)
  ####3 check datasets
  if( (dim(abar)[1] == (p1 + p2) ) ) 
  {print("We are good to go!")} else 
  { print("Wrong Datasets (X1 X2 dimensions not match with ones used in SmCCNet)!") }
  # feature names 
  rna_names <- colnames(X1)
  micro_names <- colnames(X2)
  fea_names <- colnames(abar)
  ############### including the module 0 ############ 
  module_s_0 <- add_module_0(abar, modules)
  ########## correlation against the pheotype ###########
  # the rescaled value for each module
  fea_list <- levels_mo_feature(X1 , 
                                X2 , 
                                module_s_0)
  ### the stacked dataframe of pearson correlation against the CD14 ##########
  modules_Pearson <- corr_list_module_NA (fea_list, Y, 0:n_networks, modulelabel = "module")
  ########### get correlation for different edges cuts #########
  adj_cut <- vector("list", length(edges))
  features_cut <- vector("list", length(edges))
  trimmed_Pearson <- vector("list", length(edges))
  for (i in 1: length(edges)){
    edge_cut <- edges[i]
    ############### get the adjacency matrix ##############
    # tmp <- signed_sim_matrix_cut(X1 = filtered_rlog[, filtered_outlier], 
    #                              X2 = mibi[, mibi_outlier], abar,
    #                              modules, edgecut = edge_cut)
    tmp <- signed_sim_matrix_cut(X1 , 
                                 X2 , abar,
                                 modules, edgecut = edge_cut)
    adj_cut[[i]] <- tmp
    ############### get feature levels  ##############
    tmp2 <- levels_mo_feature_mat(X1 ,
                                  X2 , tmp)
    features_cut[[i]] <- tmp2
    ######## pearson correlation against sCD14 ##########
    
    # trimmed_Pearson[[i]] <- corr_list_module_noNA(tmp2, Y, 1:n_networks, 
    #                          modulelabel = paste(edge_cut, "trimmed module"), n_na_LPS)
    trimmed_Pearson[[i]] <- corr_list_module(tmp2, Y, 1:n_networks, 
                                             modulelabel = paste(edge_cut, "trimmed module"))
  }
  ########## focus on edge 0.1#########
  message("cor.test of 0.1 edge cut nodes!")
  n = which(edges %in% 0.1)
  features_01 = features_cut[[n]][lapply(features_cut[[n ]],ncol)>0]
  if( (max( Y[,1], na.rm = TRUE) == 1) &  
          (min( Y[,1], na.rm = TRUE) == 0)){
    
    ######### see pca_getPC1.R #########
    
  } else{
    for (i in length(features_01)) {
      features_01_tmp = features_01[[i]] %>% as.data.frame()
      tmp_cor = cor_test_df(features_01_tmp, Y[,1]) %>% round(., 4) %>% as.data.frame() %>% 
        t() %>% as.data.frame() 
      id = sim_micro_names(rownames(tmp_cor)) 
      final01_cor <-  tmp_cor  %>% dplyr::mutate(`FDR (correlation)` = 
                                                   round( p.adjust(pvalue, method = "BH"), 4),
                                                 id = id) %>%
        plyr::arrange(pvalue) %>% dplyr::rename(`p (correlation)` = pvalue) %>%
                              dplyr::select(id, `FDR (correlation)`, everything() )
      rownames(final01_cor) <- id
      write.xlsx( final01_cor, 
                  file = paste0(dir, run, "01nodes_against_Pheno.xlsx"))
    }
    
  }
  
    ###### combine all trimmed correlation 
  trimmed_Pearson_all <- do.call("rbind", trimmed_Pearson)
  modules_cor <- rbind(modules_Pearson, trimmed_Pearson_all)
  ########## the plot ###############
  corr_modules(modules_cor, modules_cor$ind, modules_cor$values,  "", "Pearson's r",
               run, dir)
  ########### PC1s against the phenotype ###########
  ## original ones 
  pc1_modules <- do.call(  "cbind", get_pc1_listfeature(fea_list))
  
  
  n_modules <- length(n_strong_modules)
  ## trimmed ones (PC1s)
  mat <- matrix(NA, nrow( fea_list[[n_strong_modules[1]]]), 
                length(edges)*n_modules ) 
  for (i in 1: length(edges)){
    mat[ , (i*n_modules + 1 - n_modules) :( i*n_modules )  ] <- 
      get_pc1_listfeature(features_cut[[i]]) %>% unlist
    
  }
  if (n_modules == 1){
    message("Plot PC1 against phenotype")
    df_trend = data.frame(PC1 = mat[,1], Y = Y[,1])
    p = ggplot(data = df_trend, aes(x = Y, y = PC1)) +
      geom_point()+
      geom_smooth(method=lm, se=FALSE, color="black") +
      theme_bw() +
      labs( x = Y_name , y = "PC1 (features at 0.1 edge cut)") +
      theme(axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            text=element_text(family="Arial"))
    print(p)
    ggsave(filename = paste0(run, "PC1_against_pheno.tiff"),device = NULL,
           path = dir, dpi = 300, compression = "lzw" )
    
  }else if (n_modules == 2){
    
  }
  
  ## combind all PC1s 
  pc1_cor <- cbind(pc1_modules, mat)
  names_edges <- paste(edges, "Trimmed Module" )
  all_names <- rep(NA, n_modules*length(edges))
  for( i in 1:length(names_edges) ){
    all_names[(i*n_modules + 1 - n_modules) :( i*n_modules )]  <- 
                 paste(names_edges[i], n_strong_modules)
  }
  colnames(pc1_cor) <- c(paste("Module", 0:n_networks), all_names )
  ## PC1s correlation against Phenotype ######33
  df <- cor_test_df(pc1_cor, Y[,1]) %>% round(., 3) %>% as.data.frame()
  # p-value
  df1 <- melt(df[1,], na.rm = TRUE) %>% as.data.frame() %>% set_colnames(c("Modules",
                                                                           "p-value"))
  # rearrange the data
  data1 <- df1[c( (n_networks + 1):1, 
                  (n_networks + 1 + length(edges)*n_modules): (n_networks + 1 + 1)  ) , ]
  # Pearson's r 
  df2 <- melt(df[2,], na.rm = TRUE) %>% as.data.frame() %>% set_colnames(c("Modules",
                                                                           "Pearson's r"))
  data2 <- df2[c( (n_networks + 1):1, 
                  (n_networks + 1 + length(edges)*n_modules): (n_networks + 1 + 1)  ) , ]
  pc1_sum <- cbind(data1, data2[,2], paste0(data2[,2]," (", data1[,2],  ")")) %>%
    set_colnames(c("Modules", "p value", "Pearson's r", "Pearson's  r (p value)"))
  #### save the summary ####333
  write.xlsx(pc1_sum, 
             file = paste0(dir, run, "PC1_against_Pheno.xlsx"))
  
}

############## correlation of the 0.1 strong module against phenotypes #########
pc1_other_trait <- function(abar, modules, X1, X2, run, trait_df,
                            dir, edge_cut){
  message("Please load the abar, modules and Ws first, this one needs the Y!")
  message(run)
  message("Check the dimension of microbiome dataset: ")
  p1 <-  ncol(X1); p2 <-  ncol(X2); X2 %>% ncol() %>% print()
  # number of feature
  ####3 check datasets
  if( (dim(abar)[1] == (p1 + p2) ) ) 
  {message("We are good to go!")} else 
  { message("Wrong Datasets (X1 X2 dimensions not match with ones used in SmCCNet)!") }
  ########## correlation against the pheotype ###########
  # adj_cut <- vector("list", length(edges))
  # features_cut <- vector("list", length(edges))
  ########3 the adjacency matrix ###########
  tmp <- signed_sim_matrix_cut(X1 , 
                               X2 , abar,
                               modules, edgecut = edge_cut)
  ############### get feature levels  ##############
  tmp2 <- levels_mo_feature_mat(X1 = filtered_rlog[, filtered_outlier],
                                X2 = mibi[, mibi_outlier], tmp)
  ###### check the dimension of the feature #######
  # tmp2 %>% 
  #   ###### non-zero list #######
  # keep(., negate(is_empty)) %>% as.data.frame() %>% dim()
  pc1 = get_pc1_listfeature(tmp2) %>% unlist
  ######## correlation test against other traits######
  corre_res = cor_test_df(trait_df, pc1)
  row.names(corre_res) = c("p-value", "Pearson-r")
  corre_res %<>% t()
  write.xlsx(corre_res,  
             file = paste0(dir, run, "cut", edge_cut, "_PC1_other_pheno.xlsx"), 
             sheetName = run, row.names = TRUE) 
  
  ############ generate plots ###############
  message("Plot PC1 against phenotype!")
  for(i in 1:ncol(trait_df)){
    trait = trait_df[, i]
      Y_name = colnames(trait_df)[i] 
      df_trend = data.frame(PC1 = pc1, Y = trait)
      p = ggplot(data = df_trend, aes(x = PC1, y = Y)) +
        geom_point()+
        geom_smooth(method=lm, se=FALSE, color="black") +
        theme_bw() +
        labs( x = paste("PC1 (features with 0.1 edge trimming) of case", run) , y = Y_name ) +
        theme(axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 16),
              axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              text=element_text(family="Arial"))
      print(p)
      ggsave(filename = paste0(run, "PC1_vs_", Y_name, ".tiff"),device = NULL,
             path = dir, dpi = 300, compression = "lzw" )
      
  }
  
  
  return(corre_res)
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
# corr_modules(test_modules, test_modules$ind, test_modules$values,  "", "Pearson's r",
#              run, dir)
