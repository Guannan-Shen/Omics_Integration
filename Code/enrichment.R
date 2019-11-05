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
top_hubs_bymean <- function(abar, modules_i){
  # get modules[[i]] from abar
  df = abar[modules_i, modules_i ] %>% as.matrix() %>% as.data.frame()
  # edges strength by mean 
  mean_edge = mean(apply(df, 2, mean))
  edges_nodes =  apply(df > mean_edge, 1, sum)  %>%
    as.data.frame() %>% tibble::rownames_to_column("Nodes")
  colnames(edges_nodes)[2] = "Num. Edges"
  #### save top 25% and top50%
  n25 = ceiling(0.25* nrow(edges_nodes))
  n50 = ceiling(0.5* nrow(edges_nodes))
  nodes25 = edges_nodes %>% dplyr::top_n(n25) %>% plyr::arrange(plyr::desc(`Num. Edges`) )
  nodes50 = edges_nodes %>% dplyr::top_n(n50) %>% plyr::arrange(plyr::desc(`Num. Edges`) )
  return(list(nodes25 = nodes25, nodes50 = nodes50))
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
           pvalueCutoff  = 0.05,
           qvalueCutoff  = fdrcut
           #, readable      = TRUE  # map gene_ID ensembl ID to symbol
         )
  return(overGO)
}

# kegg pathway clusterprofiler
clusterPro_kegg <- function(symbols, fdrcut){
  # gene ID transfer
  eg = bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  # using ncbi id 
  gg_test = enrichKEGG(gene         = eg$ENTREZID,
                        organism     = 'hsa',
                        keyType       = 'ncbi-geneid',
                        pAdjustMethod = "BH",
                        qvalueCutoff = fdrcut,
                        pvalueCutoff = 0.2
                       #, pvalueCutoff = 0.05
                       )
  return(gg_test)
}


######## summary of clusterprofiler ########
summary_clusterPro <- function(abar, modules, p1, ontol, fdrcut){
  n_networks = length(modules)
  features = colnames(abar)
  go_list = vector("list", n_networks)
  kg_list = vector("list", n_networks)
  # genelists 
  for ( i in 1:n_networks){
    # gene index 
    index = modules[[i]] [modules[[i]] <= p1]
    genes = features[index]
    print(i)
    print(ontol)
    go_test = clusterPro_GO( genes, ontol, fdrcut)
    go_list[[i]] =  go_test
    print(barplot(go_test, showCategory = length(go_test$qvalue) ))
    # 
    print("KEGG")
    gg_test = clusterPro_kegg(genes, fdrcut)
    kg_list[[i]] =  gg_test
    print(barplot(gg_test, showCategory = length(gg_test$qvalue) ))

    # kegg_enrichr( genes) %>% kable %>% print()
  }
  return(list(GO = go_list, KEGG = kg_list))
}
## test run #########
## summary_clusterPro(abar, modules[[5]], p1, 'CC', fdrcut = 0.2)

## for one module ######3
# index = modules[[5]] [modules[[5]] <= p1]
# genes =colnames(abar)[index]
# go_test = clusterPro_GO( genes, ontol = "CC", fdrcut = 0.2)
# gg_test = clusterPro_kegg(genes, fdrcut = 0.2)
# print(barplot(gg_test, showCategory = length(gg_test$qvalue) ))
# go_list[[i]] =  go_test
# print(barplot(go_test, showCategory = length(go_test$qvalue) ))

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

summary_networks <- function(abar, modules, p1, n){
  n_networks <- length(n)
  features <- colnames(abar)
  
  # for ( i in 1:n_networks){
  #   print("Total nodes:")
  #   print(length(modules[[i]]))
  #   top_hubs_bymean(abar, modules[[i]], n = 20) %>% kable %>% print()
  # }
  # 
  # genelists 
  df_list = vector("list", n_networks)
  for ( i in n){
    # gene index 
    index = modules[[i]] [modules[[i]] <= p1]
    genes = features[index]
    print(i)
    print("KEGG")
    kegg_enrichr( genes) %>% kable %>% print()
    kegg = kegg_enrichr( genes)
    print("GO Cellular Component")
    GOCC_enrichr( genes) %>% kable %>% print()
    GOCC = GOCC_enrichr( genes)
    df_list[[i]] = list(KEGG = kegg, GOCC = GOCC)
  }
  return(df_list)
}

###### enrichment analysis for this strong modules 
robust_module_enrich_nodes <- function(abar, modules, X1, X2, Y, run, n_strong_modules, dir){
  print("Please load the abar, modules and Ws first, this one focuses on the Strong/Robust module!")
  print(run)
  print("Check the dimension of microbiome dataset: ")
  mibi[, mibi_outlier] %>% ncol() %>% print()
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
  ###### get the gene symbol ##########
  index = modules[[n_strong_modules]] [modules[[n_strong_modules]] <= p1]
  genes = colnames(abar)[index]
  ## go and kegg  barplot ###3
  #### cluster profiler ####3
  go_test = clusterPro_GO( genes, ontol = "CC", fdrcut = 0.2)
  gg_test = clusterPro_kegg(genes, fdrcut = 0.2)
  # gene ID transfer
  # get the ENTREZ ID 
  eg = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  #### save the genelists 
  write.xlsx(eg, 
             file = paste0(dir, run, "Module",  n_strong_modules, "genelists.xlsx"))
  
  barplot(go_test, showCategory = length(go_test$qvalue) ) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          text=element_text(family="Arial"))
  
  ggsave(filename = paste0(run, n_strong_modules,"CC.tiff"),device = NULL,
         path = dir , dpi = 300, compression = "lzw" , 
         width = 10, height = 8, units = "in")
  
  barplot(gg_test, showCategory = length(gg_test$qvalue) ) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          text=element_text(family="Arial"))
  
  ggsave(filename = paste0(run, n_strong_modules,"kegg.tiff"),device = NULL,
         path = dir , dpi = 300, compression = "lzw", 
         width = 10, height = 8, units = "in" )
  #############save top nodes ############
  ####### summary of nodes and edges ##########
  #### top 25% and top 50%
  topnodes <-  top_hubs_bymean(abar, modules[[n_strong_modules]]) 
  topnodes$nodes25 %>%  
    dplyr::mutate(Nodes = sim_micro_names(Nodes)) %>% 
    write.xlsx(., file = paste0(dir, n_strong_modules,"top_0.25_nodes_strongmodules.xlsx"))
  topnodes$nodes50 %>%  
    dplyr::mutate(Nodes = sim_micro_names(Nodes)) %>% 
    write.xlsx(., file = paste0(dir, n_strong_modules,"top_0.5_nodes_strongmodules.xlsx"))
  ###### enrichr ###############
  enrich_r <- summary_networks(abar, modules, p1, n_strong_modules)
  ##3 only save the nonZero ones 
  write.xlsx(enrich_r[[n_strong_modules]], 
             file = paste0(dir, run, n_strong_modules, "enrichr_kegg_CC.xlsx"))
  
}

######### comparing nodes #############

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


