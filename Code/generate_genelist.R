options(stringsAsFactors = F)
options(dplyr.width = Inf)

library(readxl)
library(tidyverse)
library(magrittr)
setwd("~/Documents/gitlab/Omics_Integration/DataRaw/hiv_infected_un/")
getwd()

get_rlog <- function(){
  dir = "~/Documents/gitlab/Omics_Integration/DataRaw/hiv_infected_un/"
  
  ###### preprocessing of data ########
  small_lib = c("MIHIV998")
  # shared sample size across datasets
  shared_sam = c( "MIHIV124", "MIHIV132", "MIHIV138", "MIHIV154", "MIHIV178", "MIHIV255", 
                  "MIHIV278", "MIHIV286", "MIHIV323", "MIHIV361", "MIHIV391", "MIHIV404", 
                  "MIHIV428", "MIHIV493", "MIHIV582", "MIHIV594", "MIHIV648", "MIHIV683", 
                  "MIHIV708", "MIHIV716", "MIHIV819", "MIHIV825", "MIHIV839", "MIHIV914", 
                  "MIHIV947", "MIHIV972", "MIHIV998")
  # final_sam = shared_sam[ shared_sam %nin% small_lib]
  # use all 27 samples
  final_sam = shared_sam
  
  data = read.xlsx( paste(dir, "rlog_counts_linear_regression.xlsx",sep = "") ) %>% as.data.frame()
  # rename
  colnames(data)[-c(1,2)] = c(gsub("C", "MIHIV",  colnames(data)[3:15]  ),
                              gsub( "H", "MIHIV", colnames(data)[16:34] ) )
  rownames(data) <- NULL
  # check sample 
  # all the final samples are in the transcriptome data. 
  print("Check Samples, Match: ")
  print(sum(colnames((data)[-c(1,2)]) %in% final_sam) == length(final_sam))
  # subset of libraries
  # subset of samples, and observation row variable col
  df = data %>% dplyr::select( c("Gene_ID", final_sam) ) %>% column_to_rownames("Gene_ID") %>%
    as.matrix() %>% t %>% as.data.frame
  # IDs of genes
  gene_ID_Sym = data %>% dplyr::select( c("Gene_ID", "Symbol") )
  # 
  return(list(df = df, sym = gene_ID_Sym))
}


get_tmm <- function(){
  dir = "~/Documents/gitlab/Omics_Integration/DataRaw/hiv_infected_un/"
  
  ###### preprocessing of data ########
  small_lib = c("MIHIV998")
  # shared sample size across datasets
  shared_sam = c( "MIHIV124", "MIHIV132", "MIHIV138", "MIHIV154", "MIHIV178", "MIHIV255", 
                  "MIHIV278", "MIHIV286", "MIHIV323", "MIHIV361", "MIHIV391", "MIHIV404", 
                  "MIHIV428", "MIHIV493", "MIHIV582", "MIHIV594", "MIHIV648", "MIHIV683", 
                  "MIHIV708", "MIHIV716", "MIHIV819", "MIHIV825", "MIHIV839", "MIHIV914", 
                  "MIHIV947", "MIHIV972", "MIHIV998")
  # final_sam = shared_sam[ shared_sam %nin% small_lib]
  # use all 27 samples
  final_sam = shared_sam
  
  data = read.xlsx( paste(dir, "TMM_normalized_counts.xlsx",sep = "") ) %>% as.data.frame()
  # rename
  colnames(data)[-c(1,2)] = c(gsub("C", "MIHIV",  colnames(data)[3:15]  ),
                              gsub( "H", "MIHIV", colnames(data)[16:34] ) )
  rownames(data) <- NULL
  # check sample 
  # all the final samples are in the transcriptome data. 
  print("Check Samples, Match: ")
  print(sum(colnames((data)[-c(1,2)]) %in% final_sam) == length(final_sam))
  # subset of libraries
  # subset of samples, and observation row variable col
  df = data %>% dplyr::select( c("Gene_ID", final_sam) ) %>% column_to_rownames("Gene_ID") %>%
    as.matrix() %>% t %>% as.data.frame
  # IDs of genes
  gene_ID_Sym = data %>% dplyr::select( c("Gene_ID", "Symbol") )
  # 
  return(list(df = df, sym = gene_ID_Sym))
}

# #variance.

# var_tmm <- get_tmm()$df %>% apply(., 2, var) %>% as.data.frame()
# colnames(var_tmm) <- "values"
# var_tmm <- var_tmm %>% filter(values < 200)
# 
# var_rlog <- get_rlog()$df %>% apply(., 2, var) %>% as.data.frame()
# colnames(var_rlog) <- "values"
# # mean
# mean_rlog <- get_rlog()$df %>% apply(., 2, mean) %>% as.data.frame()
# colnames(mean_rlog) <- "values"
# 
# mean_tmm <- get_tmm()$df %>% apply(., 2, mean) %>% as.data.frame()
# colnames(mean_tmm) <- "values"
# mean_tmm <- mean_tmm %>% filter(values < 300)

# var_tmm <- get_tmm()$df %>% apply(., 2, var) %>% as.data.frame()
# colnames(var_tmm) <- "values"
# mean_tmm <- get_tmm()$df %>% apply(., 2, mean) %>% as.data.frame()
# colnames(mean_tmm) <- "values"

## density plot
density_values <- function(data, xlab){
  p = ggplot(data, aes(x = values)) + 
    # alpha controls the transparency 
    geom_density(alpha = 0.6) +
    # gray scale color
    # grey scale plot
    scale_colour_grey() +
    # start = 0.2, end = 0.8
    scale_fill_grey() +
    theme_bw() +
    theme(legend.position="bottom", legend.box = "horizontal") +
    labs(caption = 
           # unlist(strsplit( as.character( substitute(data) ), "_", fixed = TRUE))[2]
           as.character( substitute(data)) ,
         y = "Density",
         x = xlab) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          text=element_text(family="Arial"))
  print(p)
}

# # density_values(var_rlog)
# density_values(var_tmm, xlab = "Variances of Genes (TMM Normalized)" )
# # # density_values(mean_rlog)
# # density_values(mean_tmm, xlab = "Means of Genes (TMM Normalized)" )
# ggsave(filename =  paste0( "Variance < 200 trans", ".tiff"), device = NULL,
#        path = "~/Documents/gitlab/Omics_Integration/DataProcessed/plots/non_collapse/",
#        dpi = 300, compression = "lzw", 
#        width = 7, height = 5, units = "in")


elbow_var_tmm <- function(vars_cutoffs){
  # get the data, tmm
  data = get_tmm()$df
  # var cutoff
  cuts = vars_cutoffs 
  n_features = c(NULL)
  for(i in 1:length(cuts) ){
    n_features[i] <- data %>% select_if(~var(.) > cuts[i]) %>% ncol(.)
  }
  return(n_features )
}

elbow_mean_tmm <- function(means_cutoffs){
  # get the data, tmm
  data = get_tmm()$df
  # var cutoff
  cuts = means_cutoffs 
  n_features = c(NULL)
  for(i in 1:length(cuts) ){
    n_features[i] <- data %>% select_if(~mean(.) > cuts[i]) %>% ncol(.)
  }
  return(n_features )
}

##### elbow plot of variance cutoff
# vars_cutoffs <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
# n_features <- elbow_var_tmm(vars_cutoffs )
# coords <- paste("(", vars_cutoffs, ", ", n_features, ")", sep="")
# coords
# 
# p = ggplot(mapping =  aes(x = vars_cutoffs, y = n_features )) + 
#   theme_bw() +
#   geom_line() +
#   labs(x = paste("Variance Cutoffs"),
#        y = paste0("Number of Genes") ) + 
#   geom_label(aes(x = vars_cutoffs, y = n_features , label=coords))
# 
# print(p)
# 
# 
# means_cutoffs <- c(0,  1,  2,  3,  4,  5, 6, 7, 8, 9, 10)
# n_features <- elbow_mean_tmm(means_cutoffs )
# coords <- paste("(", means_cutoffs, ", ", n_features, ")", sep="")
# coords
# 
# p = ggplot(mapping =  aes(x = means_cutoffs, y = n_features )) + 
#   theme_bw() +
#   geom_line() +
#   labs(x = paste("Mean Cutoffs"),
#        y = paste0("Number of Genes") ) + 
#   geom_label(aes(x = means_cutoffs, y = n_features , label=coords))
# 
# print(p)
#############3 variance 2 mean 5 of TMM ###############
## 
# mean_cut <- 5
# var_cut <- 2
# filtered_trans <- get_tmm()$df %>% select_if(~mean(.) > mean_cut) %>% select_if(~var(.) > var_cut) %>% 
#             t %>% as.data.frame() %>% rownames_to_column("Gene_ID")
# dim(filtered_trans)
# # write.csv(filtered_trans, "~/Documents/gitlab/Omics_Integration/DataRaw/hiv_infected_un/filtered_transcriptome.csv", 
# #           row.names = FALSE)
# 
# filtered_rna <-  rescaled_rna(filtered_trans, rlog = T)
# filtered_rlog <- filtered_rna[[1]]
# dim(filtered_rlog)
# 
# mean_cut <- 10
# var_cut <- 5
# filtered_trans <- get_tmm()$df %>% select_if(~mean(.) > mean_cut) %>% select_if(~var(.) > var_cut) %>% 
#   t %>% as.data.frame() %>% rownames_to_column("Gene_ID")
# dim(filtered_trans)

filter_rescale_rna <- function(mean_cut, var_cut, rlog){
  filtered_trans = get_tmm()$df %>% select_if(~mean(.) > mean_cut) %>% select_if(~var(.) > var_cut) %>% 
    t %>% as.data.frame() %>% rownames_to_column("Gene_ID") %>% rescaled_rna(., rlog = rlog)
  return(filtered_trans)
}

# filtered_rlog <- filter_rescale_rna(10, 5, T)