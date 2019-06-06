## clean the transcriptome Data 
## From a TMM normalized counts and a gene list, and a list of samples 
'%nin%' <- Negate('%in%')
options(stringsAsFactors = F)

library(readxl)
library(tidyverse)
library(magrittr)
setwd("~/Documents/gitlab/Omics_Integration/DataRaw/hiv_infected_un/")
getwd()

# ###### the raw counts  ##########3
# cnts.raw <- read.delim("All_Sample_geneCounts_raw_counts.txt", header = TRUE, sep = "\t") 
# # filter
# cnts_fsym <- cnts.raw[rowSums(cnts.raw[, 4:35])>=(5*ncol(cnts.raw[, 4:35])), ]
# ## should end up around 15 - 20K genes 
# ngenes <- nrow(cnts_fsym)
# paste("The number of remaining genes: ", ngenes, sep = '')
# colnames(cnts_fsym)[-c(1:3)] <- c(gsub("C", "MIHIV",  colnames(cnts_fsym)[4:16]  ),
#                                gsub( "H", "MIHIV", colnames(cnts_fsym)[17:35] ) )
# ## from dim() we know there are 32 samples
# cnts_f <- as.matrix(cnts_fsym[, 4:35])
# rownames(cnts_f) <- cnts_fsym$Gene_ID
# 
# pheno <- data.frame(pid = colnames(cnts_f), txt = as.factor(c(rep("Control", 13), 
#                                                      rep("HIV", 19) )) )
# pheno$txt %<>% relevel("Control")
# 
# # using the function from EDASeq
# set <- newSeqExpressionSet(as.matrix(cnts_f),
#                            phenoData = data.frame(condition=as.factor(pheno$txt),
#                                                   row.names=colnames(cnts_f)))
# dds <- DESeqDataSetFromMatrix(countData = counts(set), colData = pData(set),design = ~ condition)
# dds <- estimateSizeFactors(dds)
# # dds$sizeFactor
# ## deseq normalization and DE analysis
# register(MulticoreParam(6))
# dds <- DESeq(dds)
# ##  If many of genes have large differences in counts due to the experimental design, (Based on the DE analysis)
# ## it is important to set blind=FALSE for downstream analysis.
# ## The more the size factors differ, the more residual dependence of the variance on the mean will be found in the transformed data. 
# ## rlog is a transformation which can perform better in these cases.
# rld <- rlog(dds, blind = FALSE, fitType = "parametric")
# cnts_rld <- assay(rld)


## function to get re-scaled transcriptome data for given subset of genes or whole genes
rescaled_rna <- function(genelist, rlog){
  ### rlog is whether to rlog transform or not, boolean ###
  ############ the whole genes and samples are pre-defined #########
  dir = "~/Documents/gitlab/Omics_Integration/DataRaw/hiv_infected_un/"
  ## import genelists ##
  # as.data.frame(read.delim( paste(dir, "coreISG", sep = "" )) )
  # genelist contains Gene_ID, Symbol
  ###### preprocessing of data ########
  small_lib = c("MIHIV998")
  # shared sample size across datasets
  shared_sam = c( "MIHIV124", "MIHIV132", "MIHIV138", "MIHIV154", "MIHIV178", "MIHIV255", 
                  "MIHIV278", "MIHIV286", "MIHIV323", "MIHIV361", "MIHIV391", "MIHIV404", 
                  "MIHIV428", "MIHIV493", "MIHIV582", "MIHIV594", "MIHIV648", "MIHIV683", 
                  "MIHIV708", "MIHIV716", "MIHIV819", "MIHIV825", "MIHIV839", "MIHIV914", 
                  "MIHIV947", "MIHIV972", "MIHIV998")
  final_sam = shared_sam[ shared_sam %nin% small_lib]
  # rlog or not
  if(rlog){
    print("Regularized Log Transformation will be applied!")
    # the rlog 
    data = read.xlsx( paste(dir, "rlog_counts_linear_regression.xlsx",sep = "") ) %>%
      as.data.frame() 
    
  }else{
    print("Regularized Log Transformation will Not be applied!")
    # the TMM counts
    data = read.xlsx( paste(dir, "TMM_normalized_counts.xlsx",sep = "") ) %>%
      as.data.frame()
  }
  #### subest of samples ########
  # rename
  colnames(data)[-c(1,2)] = c(gsub("C", "MIHIV",  colnames(data)[3:15]  ),
                              gsub( "H", "MIHIV", colnames(data)[16:34] ) )
  rownames(data) <- NULL
  # check sample 
  # all the final samples are in the transcriptome data. 
  print("Check Samples, Match: ")
  print(sum(colnames((data)[-c(1,2)]) %in% final_sam) == length(final_sam))
  
  # IDs of genes
  gene_ID_Sym = data %>% dplyr::select( c("Gene_ID", "Symbol") )
  # subset of samples, TMM
  df = data %>% dplyr::select( c("Gene_ID", final_sam) ) %>% column_to_rownames("Gene_ID") %>%
    as.matrix() %>% t %>% as.data.frame
  # samples subjects
  ID = base::rownames(df)
  # rescale the data to mean 0 and variance 1 at the gene level, (for each column, within Column (feature))
  if(missing(genelist)){
    print("No gene list provided, will use the whole Transcriptome")
    # rescale using the scale function
    results = data.frame(apply(df, 2, scale)) 
    rownames(results) = rownames(df)
    genes = gene_ID_Sym
  }
  else{
    # using plyr::join to make sure the order of the gene list 
    print("Use a subset of genes")
    tmp = df[ ,colnames(df) %in% genelist$Gene_ID]
    # rescale using the scale function
    results = data.frame(apply(tmp, 2, scale))
    # rownames is the subjects
    base::rownames(results) = ID
    # get list of genes symbols
    genes = gene_ID_Sym[ gene_ID_Sym$Gene_ID %in% genelist$Gene_ID, ]
  }
  # results is the rescaled data
  return(list(results, genes))
}


