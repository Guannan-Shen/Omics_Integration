######### this is to find the optimal strength for edge cut ####3
## might use elbow plot or visually #######

options(stringsAsFactors = F)
options(dplyr.width = Inf)
library(tidyverse)
library(magrittr)
library(igraph) ####### from distance matrix, similarity matrix to adjacency matrix 
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

################ cutted similarity matrix, similarity matrix to adjacency matrix ############
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

######## try edges ################# 
try_edges <- function(edges_i){
  print("Need X1, X2, abar, modules")
  n_edges = length(edges_i)
  n_nodes = rep(NA, n_edges)
  
  ######3 get n_nodes after cut ###########
  for(i in 1:n_edges){
    iedge = edges_i[i]
    # get cutted signed similarity matrix 
    new_modules = signed_sim_matrix_cut(X1 = filtered_rlog[, filtered_outlier], 
                                        X2 = mibi[, mibi_outlier], abar,
                                        modules, edgecut = iedge)
    # count all nodes for all modules in modules
    n_nodes[i] = total_nodes_matrix(new_modules)
  }
  return(n_nodes)
}

###### elbow plot ##########
elbow_edge <- function(iedges, n_nodes, folder, title){
  ## label
  coords = paste("(", iedges, ", ", n_nodes, ")", sep="")
  p = ggplot(mapping =  aes(x = iedges, y = n_nodes )) + 
    theme_bw() +
    geom_line() +
    labs(x = paste("Edge Strength Cutoff"),
         y = paste0("N Total Nodes"),
         caption = title ) + 
    geom_label(aes(iedges, n_nodes, label=coords))
  print(p)
  ggsave(filename = paste0(folder, title,
                           ".tiff"),  dpi = 300, compression = "lzw" )
}

##########  adjacency matrix to cytoscape ##########
adj_igraph <- function(trimmed_modules){
  n_nets = length(trimmed_modules)
  adjs <- vector("list", n_nets)
  for(i in 1:n_nets){}
}

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
# # https://bioconductor.org/packages/release/bioc/vignettes/RCy3/inst/doc/Overview-of-RCy3.html
# 
# createNetworkFromDataFrames(test_cytoscape$nodeData, test_cytoscape$edgeData , 
#                             title="Soluble CD14", collection="DataFrame Example")


############# test run ####################
################## edges #################
# edges_i <- c(0, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
# n_nodes <- try_edges(edges_i )
# 
# elbow_edge(edges_i, n_nodes, "~/Documents/gitlab/Omics_Integration/Reports/plots/",
#            "Soluble CD14; Global Transcriptome; Genus Level Mirobiome")
