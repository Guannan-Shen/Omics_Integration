######### this is to find the optimal strength for edge cut ####3
## might use elbow plot or visually #######

options(stringsAsFactors = F)
options(dplyr.width = Inf)
library(tidyverse)
library(magrittr)
library(igraph) 
library(RCy3)
library(WGCNA)
library(limma)
library(diffEnrich)
library(UpSetR)
####### from distance matrix, similarity matrix to adjacency matrix 
## or from distance matrix to nodes and edges

# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
source( paste0(dir, "Code/CV_caret_smcca.R") )
source( paste0(dir, "Code/put_together.R") )

#################3 begin with contour and prediction error.csv #########
######## test with 100 50 40 2 cd14 global transcriptome microbiome data 

#### show best penalty pair ##
print('default folder: "~/Documents/gitlab/Omics_Integration/DataProcessed/"')

# CVDir <- "CD14_Outlier3_Global_100_50_Genus_3_4foldCV/"
# show_top_l1l2(CVDir, 10)

################ get datasets #################
##### loading data ######
# clinical data 
dir = "~/Documents/gitlab/Omics_Integration/"
# load all datasets 
# source( paste0(dir, "Code/8_5_testing_dataset.R") )
# run smccnet
source( paste0(dir, "Code/put_together.R") )


## run SmCCNet get the data saved, run only once is enough ###
# edge Cut 0 to get full modules, using the SmCCNet wrapper ###33


# run_SmCCNet(X1 = filtered_rlog[, filtered_outlier], 
#             X2 = mibi[, mibi_outlier],
#             Y = CD14,
#             l1 = 0.2, 
#             l2 = 0.75, 
#             s1 = 0.7, 
#             s2 = 0.9, 
#             weights = NULL,
#             # n_na = n_na,
#             # NoTrait itself is to control whether to use Y or not 
#             NoTrait = FALSE,
#             bytrait = FALSE, 
#             EdgeCut = 0)


#########3 load the data generated from the above function ##########


# dir <- paste0( "~/Documents/gitlab/Omics_Integration/DataProcessed/", CVDir)
# load(paste0(dir, "SmCCNetWeights.RData"))

## then we have averaged similarity matrix abar, and full modules (a list), 
# dim(Ws) contains 1000 subsampling run 
######3 edge Cut optimize #############
## need labels in modules, need direction of correlation from corr matrix #######
total_nodes <- function(modules){
  n_networks = length(modules)
  total = 0
  for (i in 1:n_networks){
    total = total + length(modules[[i]])
  }
  return(total)
}

#####3 get nodes from matrix 
total_nodes_matrix <- function(modules){
  n_networks = length(modules)
  total = 0
  for (i in 1:n_networks){
    total = total + ncol(modules[[i]])
  }
  return(total)
}

# total_edges_mat <- function(node_edge_list){
#   n_networks = length(node_edge_list)
#   total = 0
#   for (i in 1:n_networks){
#     df = node_edge_list[[i]]
#     total = total + nrow(df$edgeData)
#   }
#   return(total)
# }

################ cutted similarity matrix, similarity matrix to adjacency matrix ############
signed_sim_matrix_cut <- function(X1, X2, abar, modules, edgecut){
  print("Might include module 0 in modules_0")
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

######## try edges ################# 
try_edges <- function(edges_i, X1, X2, modules){
  print("Need X1, X2, abar, modules")
  n_edges = length(edges_i)
  n_nodes = rep(NA, n_edges)
  
  ######3 get n_nodes after cut ###########
  for(i in 1:n_edges){
    iedge = edges_i[i]
    # get cutted signed similarity matrix 
    new_modules = signed_sim_matrix_cut(X1 , 
                                        X2 , abar,
                                        modules, edgecut = iedge)
    # count all edges for all modules in modules
    n_nodes[i] = total_nodes_matrix(new_modules)
  }
  return(n_nodes)
}

########### WGCNA from Adjacency matrix to edges and nodes dataframe for cytoscape #######
wrappper_adj_cyto <- function(adj_mat){
  wgcna_cytoscape = exportNetworkToCytoscape(
    adj_mat %>% as.matrix(), edgeFile = NULL,
    nodeFile = NULL, weighted = TRUE, threshold = 0.0,
    nodeNames = NULL, altNodeNames = NULL,  nodeAttr = NULL,
    includeColNames = TRUE)
  return(wgcna_cytoscape)
}

## return the indices of non-empty module
non_empty_n <- function(df_list){
  n = NULL
  for(i in 1:length(df_list)){
    if(  (ncol(df_list[[i]] ) != 0) & (nrow(df_list[[i]] ) != 0) ){
      n = c(n, i)
    }else {
      n = n
    }
  }
  return(as.numeric(n) )
}

# # non_empty_n (df)
# df = signed_sim_matrix_cut (filtered_rlog[, filtered_outlier],  mibi[, mibi_outlier], abar, modules, 0.1)
# node_edge = wrappper_adj_cyto(df[[n_strong_modules]])

# # edge cut 0.1, 0.2, 0.3, 0.4  # and module 5
# adj_matrix <- adj_cut[[3]][[5]]
# # ####### WGCNA #####3
# wgcna_cytoscape <- exportNetworkToCytoscape(
#   adj_matrix %>% as.matrix(), edgeFile = NULL,
#   nodeFile = NULL, weighted = TRUE, threshold = 0.0,
#   nodeNames = NULL, altNodeNames = NULL,  nodeAttr = NULL,
#   includeColNames = TRUE)
# ##3 tidy
# df_cyto <- tidy_wgcna_cyto(wgcna_cytoscape, 
#                            x1_name = rna_names, x2_name = micro_names)
# # # https://bioconductor.org/packages/release/bioc/vignettes/RCy3/inst/doc/Overview-of-RCy3.html
# #
# createNetworkFromDataFrames(df_cyto$nodes, df_cyto$edges, 
#                             title="Soluble CD14 cut 0.3", collection="sCD14")

# $nodeData  $edgeData

############## get edges n edges ##########
get_n_edges <- function(edges_i, X1, X2, modules){
  print("Need X1, X2, abar, modules")
  n_edges = length(edges_i)
  n_nodes = rep(NA, n_edges)
  
  ######3 get n_nodes after cut ###########
  for(j in 1:n_edges){
    jedge = edges_i[j]
    # cut adj mat 
    test = signed_sim_matrix_cut (filtered_rlog[, filtered_outlier],  
                                  mibi[, mibi_outlier], abar, modules, jedge)
    # get non-0 list
    n0_test = keep(test, negate(is_empty))
    total = 0
    for (i in 1: length(n0_test)){
      ####33 transfer to edge dataframe
      node_edge = wrappper_adj_cyto(n0_test[[i]])
      total = total + nrow(node_edge$edgeData)
    }
    ###
    n_nodes[j] = total
  }
  return(n_nodes)
}

###### elbow plot ##########
elbow_edge <- function(iedges, n_nodes, folder, title, type){
  ## label
  coords = paste( iedges, ", ", n_nodes, sep="")
  
  # geom_line
  p = ggplot(mapping =  aes(x = iedges, y = n_nodes )) + 
    theme_bw() +
    geom_line() +
    labs(x = paste("Edge Strength Cutoffs"),
         y = paste0("Num. Total ", type),
         caption = title ) + 
    geom_label(aes(iedges, n_nodes, label=coords), size = 4.5, label.size = 0.2) +
    scale_x_continuous(limits = c(-0.03, 0.53), 
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5) ) +
     theme_bw() +
    theme(legend.position="bottom", legend.box = "horizontal" ) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          text=element_text(family="Arial"))
  
  print(p)
  ggsave(filename = paste0(type, title, ".tiff"),device = NULL,
         path = folder , dpi = 300, compression = "lzw" )
}

################# summary of edge cut ##############
cut_get_list_corr <- function(abar, modules, X1, X2, Y, run, edges_i, dir, micro_name){
  print("Please load the abar, modules and Ws first, this is for the edge cut!")
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
  ########## edges cut elbow plots ##########3
  n_nodes <- try_edges(edges_i, filtered_rlog[, filtered_outlier],  mibi[, mibi_outlier], modules)
  elbow_edge(edges_i, n_nodes, dir,
             paste0(as.character(substitute(Y)),
                    "; Global Transcriptome;Mirobiome ", micro_name),
             "Nodes")
  ## num of edges 
  n_edges <- get_n_edges(edges_i, filtered_rlog[, filtered_outlier],  mibi[, mibi_outlier], modules)
  elbow_edge(edges_i, n_edges, dir,
             paste0(as.character(substitute(Y)),
                    "; Global Transcriptome;Mirobiome ", micro_name),
             "Edges")
}


############# get the number of the robust module ######33
robust_module <- function(abar, modules, X1, X2, Y, run, edge_cut, dir){
  print("Please load the abar, modules and Ws first, this is to get the robust module!")
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
  ###### get the adjacency matrix by this edge cut ###############
  # keep non-zero list element  :         keep(trimmed_list, negate(is_empty))
  trimmed_list <-  signed_sim_matrix_cut (filtered_rlog[, filtered_outlier], 
                                          mibi[, mibi_outlier], abar, modules, edge_cut)
  n_strong_modules <-   non_empty_n (trimmed_list )
  print(paste("Non zero modules after trimming:", n_strong_modules))
  print( paste("Total nodes of the robust module:", length(modules[[n_strong_modules]]) )  )
}

##########  adjacency matrix to cytoscape ##########
# adj_igraph <- function(trimmed_modules){
#   n_nets = length(trimmed_modules)
#   adjs <- vector("list", n_nets)
#   for(i in 1:n_nets){}
# }

# ############### export to CytoScape ###########
######### run datasets, modules loading and cutting in 9_30_2019 #########
cytoscapePing ()
cytoscapeVersionInfo ()

##  tidy WGCNA ready for cytoscape 
tidy_wgcna_cyto <- function(wgcna_cytoscape, x1_name, x2_name){
  options(stringsAsFactors = F)
  # nodes
  nodes = wgcna_cytoscape$nodeData %>% as.data.frame() %>% 
                      dplyr::rename(id = nodeName)  %>% 
                     # grouping of nodes 
                       dplyr::mutate(group = ifelse(id %in% x1_name, 
                                                    "Transcriptome",
                                                    "Microbiome")) %>%
                    ## simplify names 
                      dplyr::mutate(id = sim_micro_names(id)) %>% dplyr::select(id, group)
     
  # edges
  ## types of edges 
 df = wgcna_cytoscape$edgeData %>% as.data.frame()
 n = nrow(df)
 interaction_vec = rep(NULL, n)
 for ( i in 1:n){
   x = df[i, 1]
   y = df[i, 2]
   if( (x  %in% x1_name) & (y %in% x1_name) ){
     interaction_vec[i] = "Between Genes"
   } else if( (x  %in% x2_name) & (y %in% x2_name) ) {
     interaction_vec[i] = "Between Taxa"
   } else {
     interaction_vec[i] = "Across"
   }
 }
  # print(interaction_vec)
  ## edges 
  edges = df %>% dplyr::rename(source = fromNode, target = toNode) %>% 
                             dplyr::select(-fromAltName, -toAltName) %>% 
                 dplyr::mutate(interaction = interaction_vec) %>%
    ## simplify names 
                 dplyr::mutate(source = sim_micro_names(source),
                               target = sim_micro_names(target),
                               abs_weight = abs(weight) )
    
  # return results
  return(list(nodes = nodes, edges = edges))
}

########### run #########
# edge cut 0.1, 0.2, 0.3, 0.4  # and module 5
# adj_matrix <- adj_cut[[3]][[5]]
# # ####### WGCNA #####3
# wgcna_cytoscape <- exportNetworkToCytoscape(
#   adj_matrix %>% as.matrix(), edgeFile = NULL,
#   nodeFile = NULL, weighted = TRUE, threshold = 0.0,
#   nodeNames = NULL, altNodeNames = NULL,  nodeAttr = NULL,
#   includeColNames = TRUE)

##3 tidy
# df_cyto <- tidy_wgcna_cyto(wgcna_cytoscape, 
#                            x1_name = rna_names, x2_name = micro_names)
# # # https://bioconductor.org/packages/release/bioc/vignettes/RCy3/inst/doc/Overview-of-RCy3.html
# #
# createNetworkFromDataFrames(df_cyto$nodes, df_cyto$edges, 
#                             title="Soluble CD14 cut 0.3", collection="sCD14")
########### for cytoscape ############
make_cytoscape <- function(abar, modules, X1, X2, edge_cut, title, collection, n_strong_modules){
  # feature names 
  rna_names <- colnames(X1)
  micro_names <- colnames(X2)
  fea_names <- colnames(abar)
  # the adjacency matrix
  trimmed_list <-  signed_sim_matrix_cut (X1, 
                                          X2, abar, modules, edge_cut)
  # get the nodes and edges data frame
  node_edge <- wrappper_adj_cyto(trimmed_list[[ n_strong_modules ]])
  #  tidy the nodes and edges data 
  df_cyto <- tidy_wgcna_cyto(node_edge, 
                             x1_name = rna_names, x2_name = micro_names)
  ### save ranked edges #########
  ranked_edge <- df_cyto$edges %>% as.data.frame() %>% plyr::arrange(plyr::desc(abs_weight))
  
  across_edge <- df_cyto$edges %>% as.data.frame() %>% 
    dplyr::filter(interaction == "Across" )   %>% 
    plyr::arrange(plyr::desc(abs_weight))
  #### save #######
  write.xlsx(df_cyto$nodes %>% plyr::arrange(group, id), 
             file = paste0(dir, n_strong_modules, edge_cut, "nodes.xlsx"))
  write.xlsx(df_cyto$edges, file = paste0(dir, n_strong_modules, edge_cut, "edges.xlsx"))
  write.xlsx(ranked_edge, file = paste0(dir, n_strong_modules, edge_cut, "ranked_edges.xlsx"))
  write.xlsx(across_edge, 
             file = paste0(dir, n_strong_modules, edge_cut, "across_inter_edges.xlsx"))
  ## export to cytoscape
  createNetworkFromDataFrames(df_cyto$nodes, df_cyto$edges, 
                              title= title, collection= collection)
  
  
}
############# comparing nodes #############
nodes_two_by_two <- function(dir1, dir2, X1, X2, title1, title2, diff){
  message("X1 must be transcriptome!")
  # number of feature
  p1 <-  ncol(X1)
  p2 <-  ncol(X2)
  p = p1 + p2
  ## 
  data_a = read.xlsx(paste0(dir1, name1  ))
  data_b = read.xlsx(paste0(dir2, name2 )) 
  #### two by two of p1
  int_22mat <- function(){
    mat3 = matrix(data = 0, nrow = 2, ncol = 2) 
    colnames(mat3) <- paste( c("Yes", "No"), title1)
    rownames(mat3) <- paste( c("Yes", "No"), title2)
    return(mat3)
  }
  mat1 <- int_22mat()
  
  df1 = data_a[data_a$group == "Transcriptome", ]
  df2 = data_b[data_b$group == "Transcriptome", ]
  df = merge(df1,
             df2,by="id",all=TRUE) %>% plyr::arrange(id)
  mat1[2,2] = p1 - nrow(df)
  df = merge(df1,y = df2,by="id") %>% plyr::arrange(id)
  write.xlsx(df %>% plyr::arrange(id), 
             file = paste0(dir2, title2, diff, "overlapping_genes.xlsx"),
             row.names = TRUE  )
  
  mat1[1,1] = nrow(df)
  mat1[2,1] = nrow(df1)  - nrow(df)
  mat1[1,2] = nrow(df2)  - nrow(df)
  ######3 two by two of p2
  mat2 <- int_22mat()
  df1 = data_a[data_a$group == "Microbiome",  ]
  df2 = data_b[ data_b$group == "Microbiome", ]
  df = merge(df1,
             df2,by="id",all=TRUE)
  mat2[2,2] = p2 - nrow(df)
  
  df = merge(df1,y = df2,by="id") %>% plyr::arrange(id)
  write.xlsx(df , file = paste0(dir2, title2, diff, "overlapping_taxa.xlsx"),
             row.names = TRUE  )
  mat2[1,1] = nrow(df)
  mat2[2,1] = nrow(df1)  - nrow(df)
  mat2[1,2] = nrow(df2)  - nrow(df)
  ##### two by two of p1 + p2
  ## outer join
  mat3 <- int_22mat()
  df = merge(data_a,y = data_b,by="id",all=TRUE)
  mat3[2,2] = p - nrow(df)
  df = merge(data_a,y = data_b,by="id") %>% plyr::arrange(group.x , id)
  write.xlsx(df, file = paste0(dir2, title2, diff, "overlapping_all.xlsx"),
             row.names = TRUE  )
  mat3[1,1] = nrow(df)
  mat3[2,1] = nrow(data_a)  - nrow(df)
  mat3[1,2] = nrow(data_b)  - nrow(df)
  return(list(Trans = mat1, Micro = mat2, TwoOmics = mat3))
}

sum_fisher <- function(mat){
  test = fisher.test(mat)
  results = data.frame(p.value = test$p.value,
                       estimate = test$estimate,
                       method = test$method,
                       null.value = test$null.value,
                       alternative = test$alternative)
  return(results)
}
# # ########## igraph ####
# adj_igraph <- graph.adjacency(
#   adj_matrix %>% as.matrix(),
#   mode="undirected",
#   weighted=TRUE,
#   diag=FALSE
# )
# #
# createNetworkFromIgraph(adj_igraph,"myIgraph")



# Export the network into edge and node list files Cytoscape can read
# reduced_sim[[5]] %>% as.matrix()
# adj_igraph <- graph.adjacency(
#   reduced_sim[[5]] %>% as.matrix(),
#   mode="undirected",
#   weighted=TRUE,
#   diag=FALSE
# )
# adj_igraph[[1]]

# http://pablobarbera.com/big-data-upf/html/02a-networks-intro-visualization.html
# https://bioconductor.org/packages/release/bioc/vignettes/RCy3/inst/doc/Cytoscape-and-iGraph.html
# plot(adj_igraph)
# createNetworkFromIgraph(adj_igraph,"myIgraph")

# 
# library(WGCNA)
# test_cytoscape <- exportNetworkToCytoscape(
#   reduced_sim[[5]] %>% as.matrix(),
#   edgeFile = NULL,
#   nodeFile = NULL,
#   weighted = TRUE,
#   threshold = 0.0,
#   nodeNames = NULL,
#   altNodeNames = NULL,
#   nodeAttr = NULL,
#   includeColNames = TRUE)
# 
# nodes <- test_cytoscape$nodeData
# 
# 
# write.csv(test_cytoscape$edgeData, row.names = F, 
#           "~/Documents/gitlab/Omics_Integration/DataProcessed/networks_cyto/cd14edge.csv")
# 
# library(RCy3)
# cytoscapePing ()
# cytoscapeVersionInfo ()



############# test run ####################
################## edges #################
# edges_i <- c(0, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
# n_nodes <- try_edges(edges_i )
# 
# elbow_edge(edges_i, n_nodes, "~/Documents/gitlab/Omics_Integration/Reports/plots/",
#            "Soluble CD14; Global Transcriptome; Genus Level Mirobiome")
