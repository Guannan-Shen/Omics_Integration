# After Run SmCCNet 

# Optimize Contour plot 
library(plotly)
library(reshape2)
library(openxlsx)
library(tidyverse)
library(magrittr)
library(knitr)
## generate networks 
library(igraph) 
library(RCy3)
library(WGCNA)

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

# load all datasets 
source( paste0(dir, "Code/8_5_testing_dataset.R") )

# more mibi datasets getLs_diff_microbiome.R

############################ edge cut #######################
################# load abar, modules, Ws, the product of SmCCNet ###############
# filtered_rlog[, filtered_outlier],  mibi[, mibi_outlier], abar, modules

CVDir <- "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
load(paste0(dir, "SmCCNetWeights.RData"))

run <- 'Genus_20_1_sCD14'
print(run)
########### always need X1 X2 ############3
p1 <-  ncol(filtered_rlog[, filtered_outlier])
p2 <-  ncol(mibi[, mibi_outlier])
n_networks <- length(modules)
# feature names 
rna_names <- colnames(filtered_rlog[, filtered_outlier])
micro_names <- colnames( mibi[, mibi_outlier])

mibi[, mibi_outlier] %>% ncol() %>% print()

# check datasets
ifelse(dim(abar)[1] == (p1 + p2), "We are good to go!", "Wrong Mibi Data!")

###### edge cut and through this, find the best (strongest) sub-networks ########3
edges_i <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
###
print("Need X1, X2, abar, modules")
n_nodes <- try_edges(edges_i, filtered_rlog[, filtered_outlier],  mibi[, mibi_outlier], modules)
elbow_edge(edges_i, n_nodes, dir,
           "Soluble CD14; Global Transcriptome; Genus Level Mirobiome (20 1)",
           "Nodes")
## num of edges 
n_edges <- get_n_edges(edges_i, filtered_rlog[, filtered_outlier],  mibi[, mibi_outlier], modules)
elbow_edge(edges_i, n_edges, dir,
           "Soluble CD14; Global Transcriptome; Genus Level Mirobiome (20 1)",
           "Edges")

########### pick an edge cutoff ############
edge_cut <- 0.1 ######### 0.4 is for cytoscape 
trimmed_list <-  signed_sim_matrix_cut (filtered_rlog[, filtered_outlier], 
                                        mibi[, mibi_outlier], abar, modules, edge_cut)
n_strong_modules <-   non_empty_n (trimmed_list )

paste("Non zero modules after trimming:", n_strong_modules)
# keep non-zero list element  :         keep(trimmed_list, negate(is_empty))

#######33 from here focus on non-zero modules #########
###### enrichment analysis for this strong modules 
index = modules[[n_strong_modules]] [modules[[n_strong_modules]] <= p1]
genes = colnames(abar)[index]
## go and kegg
#### cluster profiler
go_test = clusterPro_GO( genes, ontol = "CC", fdrcut = 0.2)
gg_test = clusterPro_kegg(genes, fdrcut = 0.2)

barplot(go_test, showCategory = length(go_test$qvalue) )
ggsave(filename = paste0(run,"CC.tiff"),device = NULL,
       path = dir , dpi = 300, compression = "lzw" )
barplot(gg_test, showCategory = length(gg_test$qvalue) )
ggsave(filename = paste0(run,"kegg.tiff"),device = NULL,
       path = dir , dpi = 300, compression = "lzw" )
###### enrichr
enrich_r <- summary_networks(abar, modules, p1, n_strong_modules)
##3 only save the nonZero ones 
write.xlsx(enrich_r[[n_strong_modules]], 
           file = paste0(dir, run, "enrichr_kegg_CC.xlsx"))

############ get nodes, edges matrix for cytoscape  ##########
edge_cut <- 0.4 ######### 0.4 is for cytoscape 
trimmed_list <-  signed_sim_matrix_cut (filtered_rlog[, filtered_outlier], 
                                        mibi[, mibi_outlier], abar, modules, edge_cut)

node_edge <- wrappper_adj_cyto(trimmed_list[[ n_strong_modules ]])

######### run datasets, modules loading and cutting in 9_30_2019 #########
cytoscapePing ()
cytoscapeVersionInfo ()

df_cyto <- tidy_wgcna_cyto(node_edge, 
                           x1_name = rna_names, x2_name = micro_names)
df_cyto$edges %>% nrow()
# # https://bioconductor.org/packages/release/bioc/vignettes/RCy3/inst/doc/Overview-of-RCy3.html
#
createNetworkFromDataFrames(df_cyto$nodes, df_cyto$edges, 
                            title= "New Soluble CD14 cut 0.4", collection="sCD14")

write.xlsx(df_cyto$nodes, file = paste0(dir, edge_cut, "nodes.xlsx"))
write.xlsx(df_cyto$edges, file = paste0(dir, edge_cut, "edges.xlsx"))
####### summary of nodes and edges ##########
n = 20
top_hubs_bymean(abar, modules[[n_strong_modules]], n = n) %>%  
                dplyr::mutate(Nodes = sim_micro_names(Nodes)) %>% 
                     write.xlsx(., file = paste0(dir, n, "top_nodes_strongmodules.xlsx"))

############# correlation against the clinical phenotypes ##############


#################### better contour plots ###############

setwd("~/Documents/gitlab/Omics_Integration/DataProcessed/")
top <- "~/Documents/gitlab/Omics_Integration/DataProcessed/"
getwd()

########## LPS 20% 1% Genus ##########
CVDir <- "CD14_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
optimize_contour(CVDir, l1_max = 0.8, l2_max = 0.8)
  
CVDir <-  "_Unclassified_Genus_Global_100_50_20_1_4_4foldCV/"
optimize_contour(CVDir, l1_max = 0.8, l2_max = 0.8)

CVDir <-   "HIV_Unclassified_Genus_Global_100_50_20_1_3_4foldCV/"
optimize_contour(CVDir, l1_max = 0.8, l2_max = 0.8)

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

######### family ###########
CVDir <- "CD14_Unclassified_Family_Global_100_50_20_3_1_4foldCV/"
optimize_contour(CVDir, l1_max = 0.55, l2_max = 0.6)

CVDir <- "CD14_Unclassified_Family_Global_100_50_40_1_1_4foldCV/"
optimize_contour(CVDir, l1_max = 0.8, l2_max = 0.8)