options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()
library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(magrittr)
library(org.Hs.eg.db)
library(enrichR)
library(enrichplot)
# Carlson M (2019). org.Hs.eg.db: Genome wide annotation for Human

######## find most connected hubs #############
# transform to dataframe

top_hubs_bymean <- function(abar, modules_i, n){
  # get modules[[i]] from abar
  df = abar[modules_i, modules_i ] %>% as.matrix() %>% as.data.frame()
  # edges strength by mean 
  mean_edge = mean(apply(df, 2, mean))
  edges_nodes =  apply(df > mean_edge, 1, sum)  %>%
    as.data.frame() %>% tibble::rownames_to_column("Nodes")
  colnames(edges_nodes)[2] = "n_strong_edges"
  final = edges_nodes %>% dplyr::top_n(n) %>% plyr::arrange(plyr::desc(n_strong_edges) )
  return(final)
}


# GO over-representaion test of Cellular Component
clusterPro_GO <- function(symbols, ontol, fdrcut) {
  overGO = 
      enrichGO(gene      = symbols,
           # universe    = names(geneList), # background
           OrgDb         = org.Hs.eg.db,  # Genome wide annotation for Human
           keyType       = "SYMBOL", # or ENSEMBL
           ont           = ontol, # 	One of "MF", "BP", and "CC" subontologies. 
           # cellular component, biological process, molecular function 
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.01,
           qvalueCutoff  = fdrcut
           #, readable      = TRUE  # map gene_ID ensembl ID to symbol
         )
  return(overGO)
}

eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
enrichKEGG(gene         = gene,
           organism     = 'hsa',
           keyType       = 'ncbi-geneid',
           pvalueCutoff = 0.05)

######## summary of clusterprofiler ########
summary_clusterPro <- function(abar, modules, p1, ontol, fdrcut){
  print(dim(abar)) # similarity matrix
  n_networks = length(modules)
  features = colnames(abar)
  
  for ( i in 1:n_networks){
    print("Total nodes:")
    print(length(modules[[i]]))
    top_hubs_bymean(abar, modules[[i]], n = 20) %>% kable %>% print()
  }
  
  # genelists 
  
  for ( i in 1:n_networks){
    # gene index 
    index = modules[[i]] [modules[[i]] <= p1]
    genes = features[index]
    print(i)

    print("GO Cellular Component")
    clusterPro_GO( genes, ontol, fdrcut) %>% kable %>% print()
    
    # print("KEGG")
    # kegg_enrichr( genes) %>% kable %>% print()
  }
  
}

# edo <- enrichDGN(de)
# library(enrichplot)
# barplot(edo, showCategory=20)
p1 = ncol(filtered_rlog[, filtered_outlier])
index = modules[[5]] [modules[[5]] <= p1]
features = colnames(abar)
genes = features[index]
cc_test <- clusterPro_GO( genes, "CC", 0.2)
barplot(cc_test, showCategory = length(cc_test$qvalue) )

cc_test$qvalue

eg = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

gg_test <- enrichKEGG(gene         = eg$ENTREZID,
           organism     = 'hsa',
           keyType       = 'ncbi-geneid',
           pAdjustMethod = "BH",
           qvalueCutoff = 0.2,
           pvalueCutoff = 0.05)
gg_test
barplot(gg_test, showCategory = length(gg_test$qvalue) )

############# enrich r ############
# KEGG over-representaion test
# hsa	Homo sapiens (human)
kegg_enrichr <- function(genelist){
  #########3 CpGtable a string and database a vector ############## 
  ######## gene_col a string
  # the usage of call "$" and eval 
  pathways = enrichr( genelist, "KEGG_2019_Human")
  # data
  kegg = data.frame( pathways[["KEGG_2019_Human"]] ) %>%
    dplyr::filter(Adjusted.P.value <= 0.2) %>% dplyr::arrange(P.value) %>% 
    dplyr::select(Term, Overlap, Adjusted.P.value, Genes)
  return( kegg )
}



# dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "KEGG_2019_Human")
GOCC_enrichr <- function(genelist){
  #########3 CpGtable a string and database a vector ############## 
  ######## gene_col a string
  # the usage of call "$" and eval 
  pathways = enrichr( genelist, "GO_Cellular_Component_2018")
  # data
  kegg = data.frame( pathways[["GO_Cellular_Component_2018"]] )%>%
    dplyr::filter(Adjusted.P.value <= 0.2) %>% dplyr::arrange(P.value)  %>% 
    dplyr::select(Term, Overlap, Adjusted.P.value, Genes)
  return( kegg )
}

# kegg_enrichr(genes)
# GOCC_enrichr(genes)

##############33 test dataset ###################

# # with CD 14
# dir = "~/Documents/gitlab/Omics_Integration/DataProcessed/CD14_Outlier1_Global_100_50_Genus_1_4foldCV/"

# # without CD 14
# dir = "~/Documents/gitlab/Omics_Integration/DataProcessed/_Outlier1_Global_100_50_Genus_1_4foldCV/"

# load(paste0(dir, "SmCCNetWeights.RData"))
# dim(abar) # similarity matrix
# n_networks <- length(modules)
# features <- colnames(abar)    # gene symbol
# for (i in 1:n_networks){
#   
# }

# clusterPro_GO(features[modules[[1]]], "CC", fdrcut = 0.2)
# kegg_enrichr(features[modules[[1]]])
# GOCC_enrichr(features[modules[[1]]])



# top_hubs_bymean(abar, modules[[1]], 10)

summary_networks <- function(abar, modules, p1){
  print(dim(abar)) # similarity matrix
  n_networks <- length(modules)
  features <- colnames(abar)
  
  for ( i in 1:n_networks){
    print("Total nodes:")
    print(length(modules[[i]]))
    top_hubs_bymean(abar, modules[[i]], n = 20) %>% kable %>% print()
  }
  
  # genelists 
  
  for ( i in 1:n_networks){
    # gene index 
    index = modules[[i]] [modules[[i]] <= p1]
    genes = features[index]
    print(i)
    print("KEGG")
    kegg_enrichr( genes) %>% kable %>% print()
    print("GO Cellular Component")
    GOCC_enrichr( genes) %>% kable %>% print()
  }
  
}

