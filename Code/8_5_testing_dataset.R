rm(list = ls())
library(openxlsx)
library(broom)
library(tidyverse)
library(magrittr)
library(SmCCNet)

options(stringsAsFactors = F)
options(dplyr.width = Inf)
setwd("~/Documents/gitlab/Omics_Integration/")
getwd()

# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
# small functions, %nin%, %||%, isnothing
source( paste0(dir, "Code/small_fun.R") )
source( paste0(dir, "Code/put_together.R") )
########## import dataset ################
# clinical data 
source( paste0(dir, "Code/6_5_clean_clinical.R") )
# subset global rna-seq by mean, var
source( paste0(dir, "Code/generate_genelist.R") )
# clean and transform transcriptome data, subset of genes
source( paste0(dir, "Code/6_5_clean_transcriptome.R") )
# outliers 
source( paste0(dir, "Code/outliers.R") )
# source to get the load_filtered_micro_level function to get clr of RA
source( paste0(dir, "Code/5_29_Generate_filtered_Data_Microbiome.R") )

# diagnostic plots and tables
source( paste0(dir, "Code/ref_plots.R") )

###### load in ####
#### phenotype contains ID ####
clin <- rescaled_cli() 
CD14 <- clin %>% dplyr::select(CD14)
ID <- clin %>% dplyr::select(ID)
sex <- clin %>% dplyr::select(sex)
age <- clin %>% dplyr::select(age)

# disease status
HIV <- ifelse(clin$Group == "hiv", 1, 0) %>% as.data.frame()
colnames(HIV) <- "HIV"
#
LTA <- clin %>% dplyr::select(LTA)
#
LPS <- clin %>% dplyr::select(LPS)

######Inflammatory marker: ########
# no missing value 
# clin %>% colnames()
inflama <- c("CRP", "IL6", "CD4_CD45", "CD4_blood", "CD4_tissue")

CRP <- clin %>% dplyr::select(CRP)
anyNA(CRP)
IL6 <- clin[,10] %>% as.data.frame() %>% set_colnames("IL6")
anyNA(IL6)
#####
message("It turns out in meta data, 24th and 25th are the same of CD4")
CD4_CD45 <- clin[,25] %>% as.data.frame() %>% set_colnames("CD4_CD45")
anyNA(CD4_CD45)
CD4_blood <- clin[,24] %>% as.data.frame() %>% set_colnames("CD4_blood")
anyNA(CD4_blood)
CD4_tissue <- clin[,9] %>% as.data.frame() %>% set_colnames("CD4_tissue")
anyNA(CD4_tissue)

#######3 use them all together #########
inflama <- cbind(ID, CRP, IL6, CD4_CD45, CD4_blood, CD4_tissue) %>% as.data.frame()
########
n_na_LPS <- which(is.na(LPS))
## IFNb
IFNb <- clin %>% dplyr::select(IFNb)
n_na_IFNb <-  which(is.na(IFNb))

##### Transcriptome ######
rna_isgs <- as.data.frame(read.delim( paste0(dir, "DataRaw/hiv_infected_un/coreISG")) ) %>% 
  rescaled_rna(., rlog = T)
# the isgs data
isgs_rlog <- rna_isgs[[1]]
# names
colnames(isgs_rlog) <- rna_isgs[[2]]$Symbol
######## genes beta
rna_genesbeta <- as.data.frame(read.delim( paste0(dir, 
                                                  "DataRaw/hiv_infected_un/genesbeta")) ) %>% 
  rescaled_rna(., rlog = T)
genesbeta_rlog <- rna_genesbeta[[1]]
colnames(genesbeta_rlog) <- rna_genesbeta[[2]]$Symbol

isgs_outlier <- grubbs_df(isgs_rlog, 2, 10)$fdr > 0.05
print(sum(isgs_outlier))
genesbeta_outlier <- grubbs_df(genesbeta_rlog, 2, 10)$fdr > 0.05
print(sum(genesbeta_outlier))

#### global filtered #######3
mean_cut <- 100
var_cut <- 50
dim(get_tmm()$df)
filtered_rna <- filter_rescale_rna(mean_cut, var_cut, T)
filtered_rlog <- filtered_rna[[1]]
colnames(filtered_rlog) <- filtered_rna[[2]]$Symbol
print(ncol(filtered_rlog))

filtered_outlier <- grubbs_df(filtered_rlog, 2, 10)$fdr > 0.05
print(sum(filtered_outlier))
print( colnames (filtered_rlog) [!filtered_outlier])

#############33##3 different mibi ##############
prev <- 20
ra <- 1
##### "genus" or "family"
micro_level <- "genus"
# omcis_name <- "30_2_Global_100_50_Genus"

micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)

micro_clr <- micro_data[[2]] %>% as.data.frame()
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
# remove low abundance taxa
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])
# a <- colnames (mibi)
############# 20, 1 ############
prev <- 20
# ra <- 1.1
ra <- 1
micro_level <- "family"
# omcis_name <- "30_2_Global_100_50_Genus"

micro_data <- load_filtered_micro_level_samples(micro_level,  
                                                prevalence = prev, RA = ra, wd = "Ubuntu",
                                                collapse = FALSE)
micro_clr <- micro_data[[2]] %>% as.data.frame()
# colnames(micro_clr)
# rescale to mean 0 and variance 1 
mibi <- rescale_microbiome(micro_clr)
mibi <- mibi %>% dplyr::select(-Other)
mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
print(ncol(mibi))
print(sum(mibi_outlier))
print( colnames (mibi) [!mibi_outlier])
# colnames (mibi) [colnames (mibi) %nin% a]
###### summary ############
print(paste0("Transcriptome cutoffs mean: ", mean_cut, " variance: ", var_cut, 
             " Microbiome cutoffs prevalence: ", prev, " RA: ", ra, " ", micro_level ) )


#################### summary of prevalence and RA, per Taxon ################
get_prev_ra <- function(micro_data_ra){
  data = micro_data_ra
  taxa = colnames(data)
  prev = apply(data, 2, function(x){
    ( ( sum(x!=0)/nrow(data) ) %>% round(., 4) ) *100
  })
  max_ra = apply(data, 2, function(x){
    ( ( max(x) ) %>% round(., 4) ) 
  })
  df = data.frame(Taxa = taxa, Prevalence = prev, RA = max_ra) %>%  
    dplyr::filter(Taxa !=  "Other")  %>%    dplyr::mutate(
      `60% 3%` = ifelse ( (Prevalence >= 60) & (RA >= 3), "60% 3%", "N"), 
      `60% 1%` = ifelse ( (Prevalence >= 60) & (RA >= 1), "60% 1%", "N"),
      `40% 3%` = ifelse ( (Prevalence >= 40) & (RA >= 3), "40% 3%", "N"),
      `40% 1%` = ifelse ( (Prevalence >= 40) & (RA >= 1), "40% 1%", "N"),
      `20% 3%` = ifelse ( (Prevalence >= 20) & (RA >= 3), "20% 3%", "N"),
      `20% 1%` = ifelse ( (Prevalence >= 20) & (RA >= 1), "20% 1%", "N")
    )
  return(df)
}

sim_micro_names_colon <- function(micro_names){
  p2 = length(micro_names)
  micro_names_sim = rep(NULL, p2)
  
  for (i in 1:p2){
    new_vector = unlist(strsplit(micro_names[i], fixed = T, split = ":"))
    test_name = new_vector[-c(1,2)] %>% paste(., collapse = ".")
    if (length(new_vector) == 1 ) {
      micro_names_sim[i] = micro_names[i]
    } else if ( length(new_vector) == 2 ){  
      micro_names_sim[i] <-  new_vector[-1] %>%
        paste(., collapse = ".")
    } else if(test_name == 2){
      micro_names_sim[i] <- "4C0d-2"
    }
    else{
      micro_names_sim[i] <-  new_vector[-c(1,2)] %>%
        paste(., collapse = ".")}
  }
  return(micro_names_sim)
}

# df = get_prev_ra(micro_data[[1]]) %>% dplyr::mutate(Taxa = sim_micro_names_colon(df$Taxa) )
# df %>% dplyr::filter(`40% 1%` == "N")
# write.xlsx(df %>% dplyr::filter(`40% 1%` == "N"), 
#            file = paste0(dir,  "only_20_1_genus.xlsx"))
# 
# df %>% dplyr::filter((`60% 1%` == "N") & (`40% 1%` == "40% 1%") )
# write.xlsx(df %>% dplyr::filter((`60% 1%` == "N") & (`40% 1%` == "40% 1%") ), 
#            file = paste0(dir,  "upto_40_1_genus.xlsx"))
# 
# write.xlsx(df %>% dplyr::filter( (`60% 1%` != "N" )  ), 
#            file = paste0(dir,  "upto_60_1_genus.xlsx"))
# write.xlsx(df, 
#            file = paste0(dir,  "genus_taxa.xlsx"))
# nrow(df)
# 33 taxa are shared across 60 1 and 40 3 genus
# sum( (df$`60% 1%` != "Not")  &    (df$`40% 3%` != "Not")  ), 


######################## descriptive ################
# ########## summary of n classified taxa ##########
# prev_list <- seq(10, 70, 10)
# ra_list <- c(0, 0.5, 1, 2, 3, 4, 5)
# data <- matrix(NA, nrow = length(prev_list)*length(ra_list), ncol = 3)
# 
# for (i in 1:length(ra_list) ){
#   micro_level <- "family"
#   ra = ra_list[i]
#   for (j in 1:length(prev_list) ){
#     prev = prev_list[j]
#     try(
#       micro_data <- load_filtered_micro_level_samples(micro_level,  
#                                                     prevalence = prev, RA = ra, wd = "Ubuntu",
#                                                     collapse = FALSE)
#       )
#     micro_clr <- micro_data[[2]] %>% as.data.frame()
#     mibi <- rescale_microbiome(micro_clr)
#     mibi <- mibi %>% dplyr::select(-Other)
#     mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
#     n_taxa <- sum(mibi_outlier)
#     ## assign value 
#     
#     data[(i-1)*length(ra_list) + j, ] <- c(paste0(ra, "%"), paste0(prev, "%"), n_taxa)
#   }
# }
# data %>% set_colnames(c("RA", "Prevalence", "N")) %>% as.data.frame() %>% 
#                 write.csv("~/Documents/gitlab/Omics_Integration/DataProcessed/genus_nocollapse.csv",
#                           row.names =  F)
# data %>% set_colnames(c("RA", "Prevalence", "N")) %>% as.data.frame() %>%
#               tidyr::pivot_wider(., names_from = RA, values_from = N) %>% as.data.frame()  %>%
# write.csv("~/Documents/gitlab/Omics_Integration/DataProcessed/genus_nocollapse_spread.csv",
#                                                   row.names =  F)
# 
# df <- data %>% set_colnames(c("RA", "Prevalence", "N")) %>% as.data.frame() %>%
#   tidyr::pivot_wider(., names_from = RA, values_from = N) %>% as.data.frame()
# 
# ##########
# 
# df[6:7,4]  <- c(29, 25)
# df
# micro_cutoffs_prev(10, "family", c(0, 0.5, 1, 2, 3, 4, 5), as.numeric(df[1, 2:8]))
# micro_cutoffs_prev(20, "family", c(0, 0.5, 1, 2, 3, 4, 5), as.numeric(df[2, 2:8]))
# micro_cutoffs_prev(40, "family", c(0, 0.5, 1, 2, 3, 4, 5), as.numeric(df[4, 2:8]))
# micro_cutoffs_prev(60, "family", c(0, 0.5, 1, 2, 3, 4, 5), as.numeric(df[6, 2:8]))
# 
# 
# micro_cutoffs_prev(10, "genus", c(0, 0.5, 1, 2, 3, 4, 5), as.numeric(df[1, 2:8]))
# micro_cutoffs_prev(20, "genus", c(0, 0.5, 1, 2, 3, 4, 5), as.numeric(df[2, 2:8]))
# micro_cutoffs_prev(40, "genus", c(0, 0.5, 1, 2, 3, 4, 5), as.numeric(df[4, 2:8]))
# micro_cutoffs_prev(60, "genus", c(0, 0.5, 1, 2, 3, 4, 5), as.numeric(df[6, 2:8]))
# ############### figures ##########
# length(prev_list)*length(ra_list)
# df
# prev <- 60
# taxa_level <- "genus"
# n_taxa <- as.numeric(df[6, 2:8])
# 
# micro_cutoffs_prev <- function(prev, taxa_level, n_taxa){
#   coords = paste(ra_list, "% ", n_taxa, sep="")
#   
#   p = ggplot(data = df, mapping =  aes(x = ra_list, y = n_taxa )) + 
#     theme_bw() +
#     geom_line() +
#     labs(x = paste("Relative Abundance Cutoff %","(", "While prevalence =",prev ,"%)", sep = " "),
#          y = paste0("N Taxa (", tools::toTitleCase(taxa_level), ")") ) + 
#     geom_label(aes(ra_list, n_taxa, label=coords), size = 4.5, label.size = 0.2)  +
#     theme(axis.text.x = element_text(size = 14),
#           axis.text.y = element_text(size = 14),
#           axis.title.x = element_text(size = 15),
#           axis.title.y = element_text(size = 15))
#   print(p)
#   title = paste0(prev, "_prev_Ntaxa_", taxa_level)
#   ggsave(filename =  paste0( title, ".tiff"), device = NULL,
#          path = "~/Documents/gitlab/Omics_Integration/DataProcessed/plots/non_collapse/",
#          dpi = 300, compression = "lzw")
#   
# }
# mibi[, "Other"]
# 
# mibi <- 
#   
#   mibi %>% dplyr::select(-Other) %>% dim()
################### outliers ################
# df <- micro_data[[1]][, "Prote:Epsil:Helicobacter"] %>% as.data.frame() %>% 
#                     dplyr::mutate(Taxa = " Proteobacteria:Epsilonproteobacteria:Helicobacter")  %>%
#                          set_colnames(c("RA", "Taxa"))
# dot_group_violin(df, df$Taxa, df$RA, "Taxa", "RA")
# dot_group_box(df, df$Taxa, df$RA, "Taxa", "RA")

# df <- filtered_rlog[, !filtered_outlier] %>% stack() %>% as.data.frame() 
# dot_group_box(df, df$ind, df$values, "Genes", "Expression Levels (rlog)")

# hist(mibi[,"Prote.Epsil.Helicobacter"])
# hist(micro_data[[1]][,"Prote:Epsil:Helicobacter"])
# micro_data[[1]] %>% colnames()
# hist(micro_data[[1]][,"Bacte:Bacte:Alistipes"])

######## hints of weights genesbeta ##########
# df5 = get_corr_1(filtered_rlog[, filtered_outlier], "Within-Transcriptome")
# summary(df5$values)
# df5[df5$values == 1, ]
# df4 = get_corr_1(mibi[, mibi_outlier], "Within-Microbiome")
# summary(df4$values)
# df4[df4$values == 1, ]

# compare_corr(mibi[, mibi_outlier], filtered_rlog[, filtered_outlier],
#              CD14, "Genus", "sCD14")
# compare_corr(mibi[, mibi_outlier], filtered_rlog[, filtered_outlier],
#              LTA, "Genus", "LTA")
# compare_corr(mibi[-n_na_LPS, mibi_outlier], 
#              filtered_rlog[-n_na_LPS, filtered_outlier],
#              LPS[-n_na_LPS, ], "Genus", "LPS")
# compare_corr(mibi[-n_na_IFNb, mibi_outlier], 
#              filtered_rlog[-n_na_IFNb, filtered_outlier],
#              IFNb[-n_na_IFNb, ], "Genus", "IFN-Beta")

# # stats::cor(clin$LPS, mibi, use = "pairwise.complete.obs", method = "pearson")

# 
# ######## isgs ############
# df1 = get_corr(mibi[-n_na, ], isgs_rlog[-n_na, ], "Microbiome:Core-ISGs")
# # stats::cor(clin$LPS, mibi, use = "pairwise.complete.obs", method = "pearson")
# df2 = get_corr(LPS[-n_na,], mibi[-n_na, ], "LPS:Microbiome")
# df3 = get_corr(LPS[-n_na,], isgs_rlog[-n_na, ], "LPS:Core-ISGs")
# # plots
# data = rbind(df1, df2, df3)

# box_values_group(data, "Pearson Correlations Summary (Core-ISGs Microbiome LPS)")
# density_values_group(data, "Pearson Correlations Summmary (Core-ISGs Microbiome LPS)")




########### select mutate ends_with ###########
###### combine unclassified to Other , at the Relative abundance level
######## check rowSums to 1, by library (subject)
##### then get the clr 

# df <- data.frame(matrix(1, ncol = 7, nrow = 4))
# colnames(df) <- c("a", "aceae", "ales", "Cyano.4C0d.2", "Bacte.Bacte.S24.7", 
#                   "Prote.Betaproteobacteria", "other")
# remove_non <- colnames(df)[4:6]
# df %>% dplyr::mutate(other = other + rowSums(df[, colnames(df)[4:6]]) + 
#                        rowSums(dplyr::select(df, ends_with("aceae")) ) +
#                        rowSums(dplyr::select(df, ends_with("ales")) ) ) %>% 
#                   dplyr::select(-remove_non, -ends_with("aceae"), -ends_with("ales"))


# G_40_2_54 <- colnames(mibi)
# 
# # using ! to not, reverse true false
# G_40_2_out <- colnames(mibi[, !mibi_outlier])
# 
# ############# or other cutoffs #########333
# prev <- 30
# ra <- 1
# # omcis_name <- "30_2_Global_100_50_Genus"
# 
# micro_data <- load_filtered_micro_level_samples("genus",  
#                                                 prevalence = prev, RA = ra, wd = "Ubuntu")
# micro_clr <- micro_data[[2]] %>% as.data.frame()
# # rescale to mean 0 and variance 1 
# mibi <- rescale_microbiome(micro_clr)
# print(ncol(mibi))
# mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
# print(sum(mibi_outlier))
# colnames(mibi)[colnames(mibi) %nin% G_40_2_54 ]
# G_30_1_66 <- colnames(mibi)
# 
# prev <- 20
# ra <- 1
# # omcis_name <- "30_2_Global_100_50_Genus"
# 
# micro_data <- load_filtered_micro_level_samples("genus",  
#                                                 prevalence = prev, RA = ra, wd = "Ubuntu")
# micro_clr <- micro_data[[2]] %>% as.data.frame()
# # rescale to mean 0 and variance 1 
# mibi <- rescale_microbiome(micro_clr)
# print(ncol(mibi))
# mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
# print(sum(mibi_outlier))
# colnames(mibi)[colnames(mibi) %nin% G_30_1_66 ]
