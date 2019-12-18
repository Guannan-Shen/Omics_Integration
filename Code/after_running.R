# After Run SmCCNet 

# Optimize Contour plot 
library(plotly)
library(knitr)
library(reshape2)
library(openxlsx)
library(tidyverse)
library(magrittr)
library(knitr)
## generate networks 
library(igraph) 
library(RCy3)
library(WGCNA)
library(diffEnrich)

options(stringsAsFactors = F)
options(dplyr.width = Inf)
# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
# small functions, %nin%, %||%, isnothing
source( paste0(dir, "Code/small_fun.R") )
# ref plots and correlation 
source( paste0(dir, "Code/ref_plots.R") )
# get module 0 and PCA
source( paste0(dir, "Code/corr_pheno.R") )
# Optimize edge cut 
source( paste0(dir, "Code/edge_cut.R") )
source( paste0(dir, "Code/put_together.R") )
source( paste0(dir, "Code/enrichment.R") )

source( paste0(dir, "Code/pca_getPC1.R") )

# load all datasets, the mibi datasets may vary
# source( paste0(dir, "Code/8_5_testing_dataset.R") )

# more mibi datasets getLs_diff_microbiome.R

############################ edge cut #######################
################# load abar, modules, Ws, the product of SmCCNet ###############
###### sCD14 genus 20 1 #############

CVDir <- "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))


###33 get the part of hierarchical clustering ######
abar_hclust(abar = abar, 290, 355)

# dd2 <- as.dist((1-abar)/2)
# # default settings of hierarchical clustering ##
# hc <- hclust(dd2, method = "complete")
# plot(hc)

##### similarity in abar from 0 to 1 ########
min(abar)
heatmap_abar_cormat(cormat = as.matrix(abar)[290:355, 290:355], 
                    TRUE, "complete", "sCD14_Genus_20_1_n2_new" )

bbar <- as.matrix(abar)*10000
bbar[bbar >= 1] <- bbar[bbar >= 1]/2000
bbar[bbar >= 1] <- rnorm(sum(bbar >= 1), 0, 0.2)
#bbar[bbar == 0] <- rnorm(sum(bbar == 0), 0, 0.2)
heatmap_abar_cormat(cormat = bbar[270:365, 270:365], 
                    TRUE, "complete", "sCD14_Genus_20_1_new" )
abar_hclust(abar = bbar, 250, 300)
# heatmap_abar_cormat(cormat = bbar, 
#                     TRUE, "complete", "sCD14_Genus_20_1_testall" )


###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
micro_name <- "(Genus 20 1)"
cut_get_list_corr(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = CD14, 
                  run = "Genus_20_1_sCD14" ,
                  edges_i, dir, micro_name )

########### always need X1 X2 ############3
########### pick an edge cutoff, and decide the robust module ############
edge_cut <- 0.1 ######### 0.4 is for cytoscape 
robust_module(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = CD14, 
                  run = "Genus_20_1_sCD14" ,
              edge_cut = 0.1, dir )
## get the correlation against the phenotype and PC1 against the phenotype ########
n_strong_modules <- 5
edges <- c(0.1, 0.2, 0.3, 0.4)
corr_against_Y(abar, modules, 
               X1 = (filtered_rlog[, filtered_outlier]), 
              X2 = mibi[, mibi_outlier], 
              Y = CD14, 
              run = "Genus_20_1_sCD14" ,
              n_strong_modules , dir, edges,
              Y_name = "sCD14 Levels")

########Correlation against other biomarkers such#########3
pc1_other_trait(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                X2 = mibi[, mibi_outlier], 
                run = "Genus_20_1_sCD14" , trait_df = inflama[, 2:6],
               dir = dir, edge_cut = 0.1) 

####### from here focus on non-zero modules #########
##### enrichment analysis #############
######33 One large figure and enrichr kegg, gene list #####
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = CD14, 
    
                       run = "Genus_20_1_sCD14" ,
                           n_strong_modules , dir, fdrcut = 0.2)

######### From here, it is time to generate the network plot #########
cytoscapePing ()
cytoscapeVersionInfo ()

############ get nodes, edges matrix for cytoscape  ##########
###### 0.4 has been done before ##########33
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.4, title = "New Soluble CD14 cut 0.4", collection="sCD14", 5)


###### two Omics genus 20 1 #############

CVDir <- "_Unclassified_Genus_Global_100_50_20_1_4_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))

run <- "Genus_20_1_Two_omics"
micro_name <- "(Genus 20 1)"
###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
## 1
cut_get_list_corr(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = NULL, 
                  run = run ,
                  edges_i, dir, micro_name)
edge_cut <- 0.1 ######### 0.3 is for cytoscape 
## 2 
robust_module(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
              X2 = mibi[, mibi_outlier], Y = NULL, 
              run = run ,
              edge_cut =edge_cut , dir )

edges <- c(0.1, 0.2, 0.3, 0.4)
corr_against_Y(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier], Y = CD14, 
               run = "Genus_20_1_TwoOmicssCD14" ,
               n_strong_modules = 4 , dir, edges)

#### get module 4 #####
n_strong_modules <- 4
## 3
## 4
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = NULL, 
                           run = run,
                           n_strong_modules , dir, fdrcut = 0.05)
######### From here, it is time to generate the network plot #########
cytoscapePing ()
cytoscapeVersionInfo ()

## 5
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.1, title = "New Two Omics cut 0.1", collection="sCD14", 4)
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.3, title = "New Two Omics cut 0.3", collection="sCD14", 4)

##############6  nodes overlapping ###############
## check X1 X2
mibi[, mibi_outlier] %>% ncol()
filtered_rlog[, filtered_outlier] %>% ncol()

X1 = (filtered_rlog[, filtered_outlier])
X2 = mibi[, mibi_outlier]

dir1 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/")

dir2 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "_Unclassified_Genus_Global_100_50_20_1_4_4foldCV/")
name1 <- "50.1nodes.xlsx"
name2 <- "40.1nodes.xlsx"

title1 <- "sCD14"
title2 <- "Omics Only"
diff <- "compare_sCD14_Omics"

######### 2 by 2 tables and Fisher Exact Test ##########
lapping_nodes <-  nodes_two_by_two(dir1, dir2, X1 = (filtered_rlog[, filtered_outlier])
                                     , X2 = mibi[, mibi_outlier]
                                     , title1, title2, diff = diff)
lapping_nodes$fisher_TwoOmics <- sum_fisher(lapping_nodes$TwoOmics)
lapping_nodes$fisher_Micro <- sum_fisher(lapping_nodes$Micro )
lapping_nodes$fisher_Trans <- sum_fisher(lapping_nodes$Trans ) 
write.xlsx(lapping_nodes, file = paste0(dir, 0.1, diff, ".xlsx"),
           row.names = TRUE  )

## differential enrichment analysis #
kegg_hsa <- get_kegg('hsa')
######## list 1, enrichment by diffEnrich ###########
list1 <- read.xlsx( paste0(dir1, "Genus_20_1_sCD14Module5genelists.xlsx") )
list1_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list1$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list1_pe$enrich_table , 
           file = paste0(dir1, "diffEnrich_kegg.xlsx"))

########### list 2, enrichment by diffEnrich #########
list2 <- read.xlsx( paste0(dir2, "Genus_20_1_Two_omicsModule4genelists.xlsx") )

print("Num of genes in total: ")
nrow(list2)

list2_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list2$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list2_pe$enrich_table , 
           file = paste0(dir2, "diffEnrich_kegg.xlsx"))
### 
print("Num of pathways in total: ")
nrow(list2_pe$enrich_table)
## finally ##
diff_enrich <- diffEnrich(list1_pe = list1_pe, list2_pe = list2_pe,
                          method = "none", cutoff = 0.05)

write.xlsx(diff_enrich$de_table , 
           file = paste0(dir2, diff, "diffEnrich.xlsx"))

plotFoldEnrichment(de_res = diff_enrich,pval = 0.05, N = 3 )
ggsave(filename = paste0(diff, ".tiff"),device = NULL,
       path = dir2 , dpi = 300, compression = "lzw", 
       width = 7, height = 5, units = "in" )

##################### HIV  20% 1% genus ###############

CVDir <- "HIV_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))

run <- "Genus_20_1_HIV"
micro_name <- "(Genus 20 1)"
Y <- HIV
###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
## 1
cut_get_list_corr(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = Y, 
                  run = run ,
                  edges_i, dir, micro_name)

edge_cut <- 0.1 ######### 0.3 is for cytoscape 
## 2 
robust_module(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
              X2 = mibi[, mibi_outlier], Y = Y, 
              run = run ,
              edge_cut =edge_cut , dir )
#### get module 1, 3 #####
n_strong_modules <- 1
n_strong_modules <- 3
## 3

edges <- c(0.1,  0.3)
corr_against_Y(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier], Y = CD14, 
               run = "Genus_20_1_HIV_sCD14" ,
               n_strong_modules = c(1,3) , dir, edges)

## 4
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           n_strong_modules , dir, fdrcut = 0.2)

robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           1, dir, fdrcut = 0.2)

######### From here, it is time to generate the network plot #########
cytoscapePing ()
cytoscapeVersionInfo ()

## 5
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.1, title = "New HIV cut 0.1 Module 3", collection="sCD14", 
               n_strong_modules = 3)
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.1, title = "New HIV cut 0.1 Module 2", collection="sCD14", 
               n_strong_modules = 1)

##############6  nodes overlapping ###############
dir1 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/")

dir2 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "HIV_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/")

name1 <- "0.1nodes.xlsx"
name2 <- "0.1nodes.xlsx"

title1 <- "sCD14"
title2 <- "HIV"
diff <- "compare_sCD14_HIV"

######### 2 by 2 tables and Fisher Exact Test ##########
lapping_nodes <-  nodes_two_by_two(dir1, dir2, X1, X2, title1, title2)
lapping_nodes$fisher_TwoOmics <- sum_fisher(lapping_nodes$TwoOmics)
lapping_nodes$fisher_Micro <- sum_fisher(lapping_nodes$Micro )
lapping_nodes$fisher_Trans <- sum_fisher(lapping_nodes$Trans ) 
write.xlsx(lapping_nodes, file = paste0(dir, 0.1, diff, ".xlsx"),
           row.names = TRUE  )

## differential enrichment analysis #
kegg_hsa <- get_kegg('hsa')
######## list 1, enrichment by diffEnrich ###########

## list 2, enrichment by diffEnrich, not enough of genes  #

### 


######################## LPS GENUS 20% 1% ###############
CVDir <- "LPS_Unclassified_Genus_Global_100_50_20_1_3_2foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))

run <- "Genus_20_1_LPS"
micro_name <- "(Genus 20 1)"
Y <- LPS
n_na_LPS
###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
## 1
cut_get_list_corr(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = Y, 
                  run = run ,
                  edges_i, dir, micro_name)

edge_cut <- 0.1 ######### 0.3 is for cytoscape 
## 2 
robust_module(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
              X2 = mibi[, mibi_outlier], Y = Y, 
              run = run ,
              edge_cut =edge_cut , dir )
#### get module 1, 3 #####
n_strong_modules <- 1

## 3
n_strong_modules <- 1
edges <- c(0.1, 0.2)
######## put na n_na_LPS ###########
corr_against_Y(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier], Y = Y, 
               run = run ,
               n_strong_modules , dir, edges, Y_name = "LPS Levels")

############## correlation against other traits #########3
pc1_other_trait(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                X2 = mibi[, mibi_outlier], 
                run = "Genus_20_1_LPS" , trait_df = inflama[, 2:6],
                dir = dir, edge_cut = 0.1) 
## 4
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           n_strong_modules , dir, fdrcut = 0.045)
######### From here, it is time to generate the network plot #########
cytoscapePing ()
cytoscapeVersionInfo ()

## 5
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.11, title = "New LPS cut 0.11 Module 2", collection="sCD14", 
               n_strong_modules = 1)


##############6  nodes overlapping ###############
dir1 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/")

dir2 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "LPS_Unclassified_Genus_Global_100_50_20_1_3_2foldCV/")

name1 <- "50.1nodes.xlsx"
name2 <- "10.1nodes.xlsx"

title1 <- "sCD14"
title2 <- "LPS"
diff <- "compare_sCD14_LPS"

######### 2 by 2 tables and Fisher Exact Test ##########
lapping_nodes <-  nodes_two_by_two(dir1, dir2, X1 = (filtered_rlog[, filtered_outlier]), 
                                   X2 =  mibi[, mibi_outlier]
                                     , title1, title2, diff)
lapping_nodes$fisher_TwoOmics <- sum_fisher(lapping_nodes$TwoOmics)
lapping_nodes$fisher_Micro <- sum_fisher(lapping_nodes$Micro )
lapping_nodes$fisher_Trans <- sum_fisher(lapping_nodes$Trans ) 
write.xlsx(lapping_nodes, file = paste0(dir, 0.1, diff, ".xlsx"),
           row.names = TRUE  )

## differential enrichment analysis #
kegg_hsa <- get_kegg('hsa')
######## list 1, enrichment by diffEnrich ###########
list1 <- read.xlsx( paste0(dir1, "Genus_20_1_sCD14Module5genelists.xlsx") )
list1_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list1$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list1_pe$enrich_table , 
           file = paste0(dir1, "diffEnrich_kegg.xlsx"))

## list 2, enrichment by diffEnrich, not enough of genes  #
list2 <- read.xlsx( paste0(dir2, "Genus_20_1_LPSModule1genelists.xlsx") )

print("Num of genes in total: ")
nrow(list2)

list2_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list2$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list2_pe$enrich_table , 
           file = paste0(dir2, "diffEnrich_kegg.xlsx"))
### 
print("Num of pathways in total: ")
nrow(list2_pe$enrich_table)
## finally ##
diff_enrich <- diffEnrich(list1_pe = list1_pe, list2_pe = list2_pe,
                          method = "none", cutoff = 0.05)
summary(diff_enrich)
write.xlsx(diff_enrich$de_table , 
           file = paste0(dir2, diff, "diffEnrich.xlsx"))

plotFoldEnrichment(de_res = diff_enrich,pval = 0.05, N = 10 )
ggsave(filename = paste0(diff, ".tiff"),device = NULL,
       path = dir2 , dpi = 300, compression = "lzw", 
       width = 7, height = 5, units = "in" )

################# Genus 20% 1 % ###########
################### LTA ##############

CVDir <- "LTA_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))

run <- "Genus_20_1_LTA"
micro_name <- "(Genus 20 1)"
Y <- LTA

###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
## 1
cut_get_list_corr(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = Y, 
                  run = run ,
                  edges_i, dir, micro_name)

edge_cut <- 0.1 ######### 0.3 is for cytoscape 
## 2 
robust_module(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
              X2 = mibi[, mibi_outlier], Y = Y, 
              run = run ,
              edge_cut =edge_cut , dir )
#### get module 1, 3 #####
n_strong_modules <- 1

## 3
n_strong_modules <- 1
edges <- c(0.1, 0.2, 0.3)
corr_against_Y(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier], Y = Y, 
               run = run ,
               n_strong_modules , dir, edges, Y_name = "LTA Levels")
############## correlation against other traits #########3
pc1_other_trait(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                X2 = mibi[, mibi_outlier], 
                run = "Genus_20_1_LTA" , trait_df = inflama[, 2:6],
                dir = dir, edge_cut = 0.1) 
## 4
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           n_strong_modules , dir, fdrcut = 0.1)
######### From here, it is time to generate the network plot #########
cytoscapePing ()
cytoscapeVersionInfo ()

## 5
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.2, title = "New LTA cut 0.2 Module 1", collection="sCD14", 
               n_strong_modules = 1)


##############6  nodes overlapping ###############
dir1 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/")

dir2 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "LTA_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/")

name1 <- "0.1nodes.xlsx"
name2 <- "10.1nodes.xlsx"

title1 <- "sCD14"
title2 <- "LTA"
diff <- "compare_sCD14_LTA"

######### 2 by 2 tables and Fisher Exact Test ##########
lapping_nodes <-  nodes_two_by_two(dir1, dir2, X1, X2, title1, title2)
lapping_nodes$fisher_TwoOmics <- sum_fisher(lapping_nodes$TwoOmics)
lapping_nodes$fisher_Micro <- sum_fisher(lapping_nodes$Micro )
lapping_nodes$fisher_Trans <- sum_fisher(lapping_nodes$Trans ) 
write.xlsx(lapping_nodes, file = paste0(dir, 0.1, diff, ".xlsx"),
           row.names = TRUE  )

## differential enrichment analysis #
kegg_hsa <- get_kegg('hsa')
######## list 1, enrichment by diffEnrich ###########
list1 <- read.xlsx( paste0(dir1, "Genus_20_1_sCD14Module5genelists.xlsx") )
list1_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list1$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list1_pe$enrich_table , 
           file = paste0(dir1, "diffEnrich_kegg.xlsx"))

## list 2, enrichment by diffEnrich, not enough of genes  #
list2 <- read.xlsx( paste0(dir2, "Genus_20_1_LTAModule1genelists.xlsx") )

print("Num of genes in total: ")
nrow(list2)

list2_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list2$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list2_pe$enrich_table , 
           file = paste0(dir2, "diffEnrich_kegg.xlsx"))
### 
print("Num of pathways in total: ")
nrow(list2_pe$enrich_table)
## finally ##
diff_enrich <- diffEnrich(list1_pe = list1_pe, list2_pe = list2_pe,
                          method = "none", cutoff = 0.05)
summary(diff_enrich)

################### Family 20% 1% sCD14 ################

CVDir <- "CD14_Unclassified_Family_Global_100_50_20_1_3_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))

run <- "Family_20_1_CD14"
micro_name <- "(Family 20 1)"
Y <- CD14

###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
## 1
cut_get_list_corr(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = Y, 
                  run = run ,
                  edges_i, dir, micro_name)

edge_cut <- 0.1 ######### 0.3 is for cytoscape 
## 2 
robust_module(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
              X2 = mibi[, mibi_outlier], Y = Y, 
              run = run ,
              edge_cut =edge_cut , dir )
#### get module 1, 3 #####
n_strong_modules <- 1

## 3
n_strong_modules <- 1
edges <- c(0.1, 0.2, 0.3, 0.4)
corr_against_Y(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier], Y = Y, 
               run = run ,
               n_strong_modules , dir, edges)
## 4
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           n_strong_modules , dir, fdrcut = 0.1)
######### From here, it is time to generate the network plot #########
cytoscapePing ()
cytoscapeVersionInfo ()

## 5
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.1, title = "Family sCD14 0.1 Module 1", collection="sCD14", 
               n_strong_modules = 1)
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.2, title = "Family sCD14 0.2 Module 1", collection="sCD14", 
               n_strong_modules = 1)


##############6  nodes overlapping ###############
dir1 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/")

dir2 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Family_Global_100_50_20_1_3_4foldCV/")

name1 <- "50.1nodes.xlsx"
name2 <- "10.1nodes.xlsx"

title1 <- "sCD14"
title2 <- "family_sCD14"
diff <- "compare_sCD14_family"

######### 2 by 2 tables and Fisher Exact Test ##########
lapping_nodes <-  nodes_two_by_two(dir1, dir2, X1, X2, title1, title2, diff)
lapping_nodes$fisher_TwoOmics <- sum_fisher(lapping_nodes$TwoOmics)
lapping_nodes$fisher_Micro <- sum_fisher(lapping_nodes$Micro )
lapping_nodes$fisher_Trans <- sum_fisher(lapping_nodes$Trans ) 
write.xlsx(lapping_nodes, file = paste0(dir, 0.1, diff, ".xlsx"),
           row.names = TRUE  )

## differential enrichment analysis #
kegg_hsa <- get_kegg('hsa')
######## list 1, enrichment by diffEnrich ###########
list1 <- read.xlsx( paste0(dir1, "Genus_20_1_sCD14Module5genelists.xlsx") )
list1_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list1$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list1_pe$enrich_table , 
           file = paste0(dir1, "diffEnrich_kegg.xlsx"))

## list 2, enrichment by diffEnrich, not enough of genes  #
list2 <- read.xlsx( paste0(dir2, "Family_20_1_CD14Module1genelists.xlsx") )

print("Num of genes in total: ")
nrow(list2)

list2_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list2$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list2_pe$enrich_table , 
           file = paste0(dir2, "diffEnrich_kegg.xlsx"))
### 
print("Num of pathways in total: ")
nrow(list2_pe$enrich_table)
## finally ##
diff_enrich <- diffEnrich(list1_pe = list1_pe, list2_pe = list2_pe,
                          method = "none", cutoff = 0.05)
summary(diff_enrich)
write.xlsx(diff_enrich$de_table , 
           file = paste0(dir2, diff, "diffEnrich.xlsx"))

plotFoldEnrichment(de_res = diff_enrich,pval = 0.05, N = 0 )
ggsave(filename = paste0(diff, ".tiff"),device = NULL,
       path = dir2 , dpi = 300, compression = "lzw", 
       width = 7, height = 5, units = "in" )


################################################### 
############### genus 40% 1% sCD14 #########
print("Genus 40% 1%, 59 Taxa")
CVDir <- "CD14_Unclassified_Genus_Global_100_50_40_1_0.34_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))

run <- "Genus_40_1_CD14"
micro_name <- "(Genus 40 1)"
Y <- CD14

###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
## 1
cut_get_list_corr(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = Y, 
                  run = run ,
                  edges_i, dir, micro_name)

edge_cut <- 0.1 ######### 0.3 is for cytoscape 
## 2 
robust_module(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
              X2 = mibi[, mibi_outlier] , Y = Y, 
              run = run ,
              edge_cut = edge_cut , dir )

#### get module 1, 2 #####
edges <- c(0.1, 0.2)
n_strong_modules <- c(1,2)
corr_against_Y(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier], Y = Y, 
               run = run ,
               n_strong_modules, dir, edges)

n_strong_modules <- 1
n_strong_modules <- 2
## 4
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           1, dir, fdrcut = 0.05)

robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           2, dir, fdrcut = 0.2)
###### make plots ###########
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.1, title = "40 1 sCD14 0.1 Module 1", collection="sCD14", 
               n_strong_modules = 1)

make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.1, title = "40 1 sCD14 0.1 Module 2", collection="sCD14", 
               n_strong_modules = 2)

##############6  nodes overlapping ###############
dir1 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/")

dir2 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_40_1_0.34_4foldCV/")
name3 <- "10.1nodes.xlsx"

name1 <- "50.1nodes.xlsx"
name2 <- "10.1nodes.xlsx"

title1 <- "sCD14"
title2 <- "40_1_sCD14"
diff <- "Module1_sCD14_prevalence"

######### 2 by 2 tables and Fisher Exact Test ##########
lapping_nodes <-  nodes_two_by_two(dir1, dir2, X1 = (filtered_rlog[, filtered_outlier])
                                     , X2 = mibi[, mibi_outlier]
                                     , title1, title2, diff)
lapping_nodes$fisher_TwoOmics <- sum_fisher(lapping_nodes$TwoOmics)
lapping_nodes$fisher_Micro <- sum_fisher(lapping_nodes$Micro )
lapping_nodes$fisher_Trans <- sum_fisher(lapping_nodes$Trans ) 
write.xlsx(lapping_nodes, file = paste0(dir, 0.1, diff, ".xlsx"),
           row.names = TRUE  )

## differential enrichment analysis #
kegg_hsa <- get_kegg('hsa')
######## list 1, enrichment by diffEnrich ###########
list1 <- read.xlsx( paste0(dir1, "Genus_20_1_sCD14Module5genelists.xlsx") )
list1_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list1$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list1_pe$enrich_table , 
           file = paste0(dir1, "diffEnrich_kegg.xlsx"))

## list 2, enrichment by diffEnrich, not enough of genes  #
list2 <- read.xlsx( paste0(dir2, "Genus_40_1_CD14Module1genelists.xlsx") )

print("Num of genes in total: ")
nrow(list2)

list2_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list2$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list2_pe$enrich_table , 
           file = paste0(dir2, diff,"diffEnrich_kegg.xlsx"))
### 
print("Num of pathways in total: ")
nrow(list2_pe$enrich_table)
## finally ##
diff_enrich <- diffEnrich(list1_pe = list1_pe, list2_pe = list2_pe,
                          method = "none", cutoff = 0.05)
summary(diff_enrich)
write.xlsx(diff_enrich$de_table , 
           file = paste0(dir2, diff, "diffEnrich.xlsx"))

plotFoldEnrichment(de_res = diff_enrich,pval = 0.05, N = 1 )
ggsave(filename = paste0(diff, ".tiff"),device = NULL,
       path = dir2 , dpi = 300, compression = "lzw", 
       width = 7, height = 5, units = "in" )


################## genus 40% 1% sCD14 ###########

CVDir <- "CD14_Unclassified_Genus_Global_100_50_40_1_1_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))

run <- "Genus_40_1_CD14"
micro_name <- "(Genus 40 1)"
Y <- CD14

###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
## 1
cut_get_list_corr(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = Y, 
                  run = run ,
                  edges_i, dir, micro_name)

edge_cut <- 0.1 ######### 0.3 is for cytoscape 
## 2 
robust_module(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
              X2 = mibi[, mibi_outlier] , Y = Y, 
              run = run ,
              edge_cut = edge_cut , dir )

#### get module 1, 3 #####
n_strong_modules <- 1
n_strong_modules <- 2
n_strong_modules <- 3

## 3
n_strong_modules <- 1
edges <- c(0.1, 0.4)
corr_against_Y(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier], Y = Y, 
               run = run ,
               n_strong_modules, dir, edges)
## 4
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           1, dir)

robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           2, dir)
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           3, dir)
######### From here, it is time to generate the network plot #########
cytoscapePing ()
cytoscapeVersionInfo ()

## 5
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.3, title = "40 1 sCD14 0.3 Module 2", collection="sCD14", 
               n_strong_modules = 2)


##############6  nodes overlapping ###############
dir1 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/")

dir2 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_40_1_1_4foldCV/")

name1 <- "0.1nodes.xlsx"
name2 <- "30.1nodes.xlsx"

title1 <- "sCD14"
title2 <- "40_1_sCD14"
diff <- "Module3_sCD14_prevalence"

######### 2 by 2 tables and Fisher Exact Test ##########
lapping_nodes <-  nodes_two_by_two(dir1, dir2, X1, X2, title1, title2)
lapping_nodes$fisher_TwoOmics <- sum_fisher(lapping_nodes$TwoOmics)
lapping_nodes$fisher_Micro <- sum_fisher(lapping_nodes$Micro )
lapping_nodes$fisher_Trans <- sum_fisher(lapping_nodes$Trans ) 
write.xlsx(lapping_nodes, file = paste0(dir, 0.1, diff, ".xlsx"),
           row.names = TRUE  )

## differential enrichment analysis #
kegg_hsa <- get_kegg('hsa')
######## list 1, enrichment by diffEnrich ###########
list1 <- read.xlsx( paste0(dir1, "Genus_20_1_sCD14Module5genelists.xlsx") )
list1_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list1$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list1_pe$enrich_table , 
           file = paste0(dir1, "diffEnrich_kegg.xlsx"))

## list 2, enrichment by diffEnrich, not enough of genes  #
list2 <- read.xlsx( paste0(dir2, "Genus_40_1_CD14Module3genelists.xlsx") )

print("Num of genes in total: ")
nrow(list2)

list2_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list2$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list2_pe$enrich_table , 
           file = paste0(dir2, diff,"diffEnrich_kegg.xlsx"))
### 
print("Num of pathways in total: ")
nrow(list2_pe$enrich_table)
## finally ##
diff_enrich <- diffEnrich(list1_pe = list1_pe, list2_pe = list2_pe,
                          method = "none", cutoff = 0.05)
summary(diff_enrich)
write.xlsx(diff_enrich$de_table , 
           file = paste0(dir2, diff, "diffEnrich.xlsx"))

plotFoldEnrichment(de_res = diff_enrich,pval = 0.05, N = 2 ) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        text=element_text(family="Arial"))


ggsave(filename = paste0(diff, ".tiff"),device = NULL,
       path = dir2 , dpi = 300, compression = "lzw", 
       width = 7, height = 5, units = "in" )

################## family 40% 1% sCD14 ###########
# list1$SYMBOL
# clusterPro_kegg_test <- function(symbols, fdrcut){
#   # gene ID transfer
#   eg = bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#   # using ncbi id 
#   gg_test = enrichKEGG(gene         = eg$ENTREZID,
#                        organism     = 'hsa',
#                        keyType       = 'ncbi-geneid',
#                        pAdjustMethod = "BH",
#                        qvalueCutoff = fdrcut,
#                        pvalueCutoff = 0.2)
#   return(gg_test)
# }
# gg_test = clusterPro_kegg_test(list1$SYMBOL, fdrcut = 0.2)
# gg_test
# barplot(gg_test, showCategory = length(gg_test$qvalue) )
# ggsave(filename = paste0(run, n_strong_modules,"kegg.tiff"),device = NULL,
#        path = dir , dpi = 300, compression = "lzw", 
#        width = 10, height = 8, units = "in" )



CVDir <- "CD14_Unclassified_Family_Global_100_50_40_1_1_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))

run <- "Genus_40_1_CD14"
micro_name <- "(Genus 40 1)"
Y <- CD14

###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
## 1
cut_get_list_corr(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = Y, 
                  run = run ,
                  edges_i, dir, micro_name)

edge_cut <- 0.1 ######### 0.3 is for cytoscape 
## 2 
robust_module(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
              X2 = mibi[, mibi_outlier] , Y = Y, 
              run = run ,
              edge_cut = edge_cut , dir )

#### get module 1, 3 #####
n_strong_modules <- 1

## 3
n_strong_modules <- 1
edges <- c(0.1, 0.4)
corr_against_Y(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier], Y = Y, 
               run = run ,
               n_strong_modules, dir, edges)
## 4
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           1, dir)


################## genus 60% 1% sCD14 #########
CVDir <- "CD14_Unclassified_Genus_Global_100_50_60_1_1_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))

run <- "Genus_60_1_CD14"
micro_name <- "(Genus 60 1)"
Y <- CD14

###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.1, 0.15,  0.25, 0.4, 0.5)
## 1
cut_get_list_corr(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                  X2 = mibi[, mibi_outlier], Y = Y, 
                  run = run ,
                  edges_i, dir, micro_name)

edge_cut <- 0.1 ######### 0.3 is for cytoscape 
## 2 
robust_module(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
              X2 = mibi[, mibi_outlier] , Y = Y, 
              run = run ,
              edge_cut = edge_cut , dir )

#### get module 1, 3 #####
n_strong_modules <- 2
n_strong_modules <- 3

## 3
edges <- c(0.1, 0.15, 0.25)
corr_against_Y(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier], Y = Y, 
               run = run ,
               n_strong_modules = c(2,3), dir, edges)
## 4
mibi[, mibi_outlier] %>% ncol()
robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           n_strong_modules = 2, dir, fdrcut = 0.2)

robust_module_enrich_nodes(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
                           X2 = mibi[, mibi_outlier],Y = Y, 
                           run = run,
                           n_strong_modules = 3, dir, fdrcut = 0.1)

######### From here, it is time to generate the network plot #########
cytoscapePing ()
cytoscapeVersionInfo ()

## 5
make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.1, title = "60 1 sCD14 0.1 Module 2", collection="sCD14", 
               n_strong_modules = 2)

make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.1, title = "60 1 sCD14 0 Module 3", collection="sCD14", 
               n_strong_modules = 3)

make_cytoscape(abar, modules, X1 = (filtered_rlog[, filtered_outlier]), 
               X2 = mibi[, mibi_outlier],
               edge_cut = 0.25, title = "60 1 sCD14 0.25 Module 2", collection="sCD14", 
               n_strong_modules = 2)


##############6  nodes overlapping ###############
dir1 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/")

dir2 <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/",
                "CD14_Unclassified_Genus_Global_100_50_60_1_1_4foldCV/")

name1 <- "50.1nodes.xlsx"
name2 <- "20.1nodes.xlsx"

title1 <- "sCD14"
title2 <- "60_1_sCD14"
diff <- "Module2_sCD14_prevalence"

######### 2 by 2 tables and Fisher Exact Test ##########
lapping_nodes <-  nodes_two_by_two(dir1, dir2, X1, X2, title1, title2, diff)


lapping_nodes$fisher_TwoOmics <- sum_fisher(lapping_nodes$TwoOmics)
lapping_nodes$fisher_Micro <- sum_fisher(lapping_nodes$Micro )
lapping_nodes$fisher_Trans <- sum_fisher(lapping_nodes$Trans ) 
write.xlsx(lapping_nodes, file = paste0(dir, 0.1, diff, "2by2.xlsx"),
           row.names = TRUE  )


## differential enrichment analysis #
kegg_hsa <- get_kegg('hsa')
######## list 1, enrichment by diffEnrich ###########
list1 <- read.xlsx( paste0(dir1, "Genus_20_1_sCD14Module5genelists.xlsx") )
list1_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list1$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list1_pe$enrich_table , 
           file = paste0(dir1, "diffEnrich_kegg.xlsx"))

## list 2, enrichment by diffEnrich, not enough of genes  #
list2 <- read.xlsx( paste0(dir2, "Genus_60_1_CD14Module3genelists.xlsx") )

print("Num of genes in total: ")
nrow(list2)

list2_pe <- pathEnrich(gk_obj = kegg_hsa, gene_list = list2$ENTREZID, cutoff = 0.2, N = 2)
write.xlsx(list2_pe$enrich_table , 
           file = paste0(dir2, diff,"diffEnrich_kegg.xlsx"))
### 
print("Num of pathways in total: ")
nrow(list2_pe$enrich_table)
## finally ##
diff_enrich <- diffEnrich(list1_pe = list1_pe, list2_pe = list2_pe,
                          method = "none", cutoff = 0.05)
summary(diff_enrich)
write.xlsx(diff_enrich$de_table , 
           file = paste0(dir2, diff, "diffEnrich.xlsx"))

plotFoldEnrichment(de_res = diff_enrich,pval = 0.05, N = 1 ) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        text=element_text(family="Arial"))


ggsave(filename = paste0(diff, ".tiff"),device = NULL,
       path = dir2 , dpi = 300, compression = "lzw", 
       width = 7, height = 5, units = "in" )


#################### better contour plots ###############

setwd("~/Documents/gitlab/Omics_Integration/DataProcessed/")
top <- "~/Documents/gitlab/Omics_Integration/DataProcessed/"
getwd()

########## LPS 20% 1% Genus ##########
CVDir <- "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
optimize_contour(CVDir, l1_max = 0.8, l2_max = 0.8)
  
CVDir <-  "_Unclassified_Genus_Global_100_50_20_1_4_4foldCV/"
optimize_contour(CVDir, l1_max = 0.5, l2_max = 0.5)

CVDir <-   "HIV_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
optimize_contour(CVDir, l1_max = 1, l2_max = 1)

CVDir <-   "IFNb_Unclassified_Genus_Global_100_50_20_1_3_3foldCV/"
optimize_contour(CVDir, l1_max = 0.8, l2_max = 0.8)

CVDir <- "LPS_Unclassified_Genus_Global_100_50_20_1_3_2foldCV/"
optimize_contour(CVDir, l1_max = 0.8, l2_max = 0.5)

CVDir <- "LTA_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
optimize_contour(CVDir, l1_max = 0.6, l2_max = 0.6)



CVDir <- "CD14_Unclassified_Family_Global_100_50_20_1_3_4foldCV/"
optimize_contour(CVDir, l1_max = 0.6, l2_max = 0.5)

CVDir <- "_Unclassified_Family_Global_100_50_20_1_4_4foldCV/"
optimize_contour(CVDir, l1_max = 0.5, l2_max = 0.5)

CVDir <- "LTA_Unclassified_Family_Global_100_50_20_1_3_4foldCV/"
optimize_contour(CVDir, l1_max = 0.8, l2_max = 0.8)

############### change cutoffs ##############

######### genus ##########33
 CVDir <- "CD14_Unclassified_Genus_Global_100_50_20_3_1_4foldCV/"
optimize_contour(CVDir, l1_max = 0.6, l2_max = 0.6)

CVDir <- "CD14_Unclassified_Genus_Global_100_50_40_1_0.34_4foldCV/"
optimize_contour(CVDir, l1_max = 0.8, l2_max = 0.8)

CVDir <- "CD14_Unclassified_Genus_Global_100_50_60_1_1_4foldCV/"
optimize_contour(CVDir, l1_max = 0.55, l2_max = 0.55)

######### family ###########
CVDir <- "CD14_Unclassified_Family_Global_100_50_20_3_1_4foldCV/"
optimize_contour(CVDir, l1_max = 0.55, l2_max = 0.6)

CVDir <- "CD14_Unclassified_Family_Global_100_50_40_1_1_4foldCV/"
optimize_contour(CVDir, l1_max = 0.8, l2_max = 0.8)