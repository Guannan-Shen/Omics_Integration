####
############# diagnostic plots ###########
## histogram and boxplot with stat_summary ##
'%nin%' <- Negate('%in%')
options(stringsAsFactors = F)
options(dplyr.width = Inf)

library(reshape2)
library(readxl)
library(openxlsx)
library(tidyverse)
library(magrittr)
library(tools)
library(wesanderson)
library(extrafont)
library(VennDiagram)

# Correlations with significance levels
# library(Hmisc)
dir = "~/Documents/gitlab/Omics_Integration/"
source( paste0(dir, "Code/put_together.R") )

setwd("~/Documents/gitlab/Omics_Integration/DataProcessed/")
getwd()


############ tidy the differential HIV untreated datasets #########
df <- read.csv("res.edger.csv") %>% dplyr::select(Gene_ID, Symbol, everything())
write.xlsx(df, "HIV_Untreated_DE_results.xlsx")

########## VennDiagram ###########
dir = "~/Documents/gitlab/Omics_Integration/DataProcessed/"
### genes ##
df = read.xlsx( paste0(dir, "side_by_side_comparison/",
                       "all_0.1_nodes_genes.xlsx"))
name = "genes"
lapping_lists <- vector('list', 4)

lapping_lists[[1]] <- list( 
             `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
              `Two Omics (NULL) (Genus 20%)` = df$TwoOmics_genus_20 %>% na.omit(),
              `HIV Status (Genus 20%)` = df$HIV_genus_20 %>% na.omit())

lapping_lists[[2]] <- list( 
              `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
              `LPS (Genus 20%)` = df$LPS_genus_20 %>% na.omit(),
              `LTA (Genus 20%)` = df$LTA_genus_20 %>% na.omit() )

lapping_lists[[3]] <- list( 
      `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
      `sCD14 (Family 20%)` = df$sCD14_family_20 %>% na.omit())

lapping_lists[[4]] <- list( 
             `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
            `sCD14 (Genus 40%)` = 
              c(as.character( df$sCD14_genus_40_1 %>% na.omit()  ),
                as.character( df$sCD14_genus_40_2 %>% na.omit()  ) ) ,
            `sCD14 (Genus 60%)` =  
              c(as.character( df$sCD14_genus_40_2 %>% na.omit()  ),
                as.character( df$sCD14_genus_40_3 %>% na.omit()  ) )  ) 



##### genes overlapping of sCD14 HIV null ####
set.seed(1) # For reproducibility of results
venn.diagram(lapping_lists[[1]], 
             filename =paste0("plots/non_collapse/", name ,"1.tiff"), 
             resolution = 300, compression = "lzw",  imagetype = "tiff",
             lwd	= rep(4,3), cex= rep(3,1), cat.cex= rep(2,1),
             cat.pos = c(180, 180, 0), 
             cat.just = list(c(0.5,-2), c(0.5,-2), c(0.5,1)) )

venn.diagram(lapping_lists[[2]], 
             filename =paste0("plots/non_collapse/", name ,"2.tiff"), 
             resolution = 300, compression = "lzw",  imagetype = "tiff",
             lwd	= rep(4,3), cex= rep(3,1), cat.cex= rep(2,1),
             cat.pos = c(180, 180, 0), 
             cat.just = list(c(0.5,-1.5), c(0.5,-1.5), c(0.5,1)) )

venn.diagram(lapping_lists[[3]], 
             filename =paste0("plots/non_collapse/", name ,"3.tiff"), 
             resolution = 300, compression = "lzw",  imagetype = "tiff",
             lwd	= rep(4,2), cex= rep(3,1), cat.cex= rep(2,1),
             cat.pos = c(0, 0), 
             cat.just = list(c(0.5,-1),  c(0.5, -1)) )

venn.diagram(lapping_lists[[4]], 
             filename =paste0("plots/non_collapse/", name ,"4.tiff"), 
             resolution = 300, compression = "lzw",  imagetype = "tiff",
             lwd	= rep(4,3), cex= rep(3,1), cat.cex= rep(2,1),
             cat.pos = c(0, 0, 180), 
             cat.just = list(c(0.6, 2), c(0.4, 2), c(0.5,0)) )

### Taxa ######
df = read.xlsx( paste0(dir, "side_by_side_comparison/",
                       "all_0.1_nodes_taxa.xlsx"))
name = "taxa"

lapping_lists[[1]] <- list( 
  `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
  `Two Omics (NULL) (Genus 20%)` = df$TwoOmics_genus_20 %>% na.omit(),
  `HIV Status (Genus 20%)` = df$HIV_genus_20 %>% na.omit())

lapping_lists[[2]] <- list( 
  `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
  `LPS (Genus 20%)` = df$LPS_genus_20 %>% na.omit(),
  `LTA (Genus 20%)` = df$LTA_genus_20 %>% na.omit() )

lapping_lists[[3]] <- list( 
  `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
  `sCD14 (Family 20%)` = df$sCD14_family_20 %>% na.omit())

lapping_lists[[4]] <- list( 
  `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
  `sCD14 (Genus 40%)` = 
    c(as.character( df$sCD14_genus_40_1 %>% na.omit()  ),
      as.character( df$sCD14_genus_40_2 %>% na.omit()  ) ) ,
  `sCD14 (Genus 60%)` =  
    c(as.character( df$sCD14_genus_40_2 %>% na.omit()  ),
      as.character( df$sCD14_genus_40_3 %>% na.omit()  ) )  ) 



##### genes overlapping of sCD14 HIV null ####
set.seed(1) # For reproducibility of results
venn.diagram(lapping_lists[[1]], 
             filename =paste0("plots/non_collapse/", name ,"1.tiff"), 
             resolution = 300, compression = "lzw",  imagetype = "tiff",
             lwd	= rep(4,3), cex= rep(3,1), cat.cex= rep(2,1),
             cat.pos = c(180, 180, 0), 
             cat.just = list(c(0.5,0), c(0.5,0), c(0.5,1.4)) )

venn.diagram(lapping_lists[[2]], 
             filename =paste0("plots/non_collapse/", name ,"2.tiff"), 
             resolution = 300, compression = "lzw",  imagetype = "tiff",
             lwd	= rep(4,3), cex= rep(3,1), cat.cex= rep(2,1),
             cat.pos = c(180, 180, 0), 
             cat.just = list(c(0.5, -0.3), c(0.5,1.5), c(0.5,1)) )

venn.diagram(lapping_lists[[3]], 
             filename =paste0("plots/non_collapse/", name ,"3.tiff"), 
             resolution = 300, compression = "lzw",  imagetype = "tiff",
             lwd	= rep(4,2), cex= rep(3,1), cat.cex= rep(2,1),
             cat.pos = c(0, 0), 
             cat.just = list(c(0.5,-1),  c(0.5, -1)) )

venn.diagram(lapping_lists[[4]], 
             filename =paste0("plots/non_collapse/", name ,"4.tiff"), 
             resolution = 300, compression = "lzw",  imagetype = "tiff",
             lwd	= rep(4,3), cex= rep(3,1), cat.cex= rep(2,1),
             cat.pos = c(0, 0, 180), 
             cat.just = list(c(0.6, 1), c(0.4, 1), c(0.5,1)) )
#### genus vs family ##
fisher.test(matrix(c(5,6,12,47), nrow = 2))

#################### individual correlation agaisnt the phenotype ########33
dim(filtered_rlog[,filtered_outlier])
dim(mibi[,mibi_outlier])
df = cbind(filtered_rlog[,filtered_outlier],
           mibi[,mibi_outlier])
dim(df)
## all features together against phenotype
cor_test_df <- function(df, y){
  mat = apply(df, 2, function(x){
    test = cor.test( x , y,  alternative = "two.sided", method = "pearson")
    res = c(pvalue = test$p.value, corr = test$estimate)
    return(res)
  })
  return(mat)
}

wilcox_test_df <- function(df, y){
  mat = apply(df, 2, function(x){
    test = wilcox.test(x, y,
                       alternative = c("two.sided"), paired = FALSE, exact = FALSE)
    diff_median = median(x[ y == 1]) - median(x[ y == 0])
    diff_mean = mean(x[ y == 1]) - mean(x[ y == 0])
    res = c(pvalue = test$p.value, 
            difference = ifelse(diff_median == 0, diff_mean, diff_median) )
    return(res)
  })
  return(mat)
}


cor_wil_test_wrapper <- function(df, Y, wilcoxon){
  if( wilcoxon){
    tmp_cor = wilcox_test_df(df, Y[,1]) %>% round(., 4) %>% as.data.frame() %>% 
      t() %>% as.data.frame() 
  }else{
    tmp_cor = cor_test_df(df, Y[,1]) %>% round(., 4) %>% as.data.frame() %>% 
      t() %>% as.data.frame() 
  }
  name = as.character(substitute(Y))
  id = sim_micro_names(rownames(tmp_cor)) 
  if( wilcoxon){
    final01_cor <-  tmp_cor  %>% dplyr::mutate(`FDR (correlation)` = 
                                                 round( p.adjust(pvalue, method = "BH"), 4),
                                               id = id) %>%
      plyr::arrange(pvalue) %>% dplyr::rename(`p (correlation)` = pvalue) %>% 
      rename(`Mean Diff.` = difference) %>%
      dplyr::select(id, `Mean Diff.`, `p (correlation)`, `FDR (correlation)`, everything() )
    
  }else {
    final01_cor <-  tmp_cor  %>% dplyr::mutate(`FDR (correlation)` = 
                                                 round( p.adjust(pvalue, method = "BH"), 4),
                                               id = id) %>%
      plyr::arrange(pvalue) %>% dplyr::rename(`p (correlation)` = pvalue) %>% 
      rename(`Pearson'r` = corr.cor) %>%
      dplyr::select(id, `Pearson'r`, `p (correlation)`, `FDR (correlation)`, everything() )
  }
  
  dir = "~/Documents/gitlab/Omics_Integration/DataProcessed/"
  write.xlsx( final01_cor, 
             file = paste0(dir, name, "_against_features.xlsx"))
  # return(final01_cor)
}

cor_wil_test_wrapper(df, CD14, wilcoxon = FALSE) 

cor_wil_test_wrapper(df, LPS, wilcoxon = FALSE)
cor_wil_test_wrapper(df, LTA, wilcoxon = FALSE)

cor_wil_test_wrapper(df, HIV, wilcoxon = TRUE)


############# read in nodes of the most concise networks ##########
strong_nodes_diff <- function(nodes, diff_file, corr_file){
  dir = "~/Documents/gitlab/Omics_Integration/DataProcessed/"
  ###### format and read the diff of genes ##########
  diff = read.xlsx( paste0(dir, diff_file) ) %>% plyr::arrange(FDR) %>% 
           dplyr::mutate(`in HIV-1 infection` = ifelse(logFC >= 0 , "Up", "Down"),
                  `p adjusted (FDR)` = format.pval(FDR, digits = 3) ) %>% 
    dplyr::select(Symbol, `in HIV-1 infection`, logFC, `p adjusted (FDR)`, PValue) 
  
  ### the nodes list #######
  df = read.xlsx( paste0(dir, 
                         nodes) )
  ## merge left join###
  df_final = base::merge(df, diff, by.x = "id",
                         by.y = "Symbol", all.x = TRUE) %>% 
                      plyr::arrange(PValue) %>% dplyr::select(-PValue)
  ############ merge data with corre against phenotype ######
  df_corr = read.xlsx( paste0(dir, corr_file) )
  df_all =  base::merge(df_corr, df_final,  by= "id") %>% 
    plyr::arrange(group, `p.(correlation)`) %>% dplyr::select(id, group, everything())
  return(df_all)
}
dir = "~/Documents/gitlab/Omics_Integration/DataProcessed/"

## sCD14 genus 20% 1% ###3
strong_nodes_diff( nodes = "side_by_side_comparison/sCD14_genus_20_1/50.4nodes.xlsx" , 
                   diff_file =  "HIV_Untreated_DE_results.xlsx",
                   corr_file = "CD14_against_features.xlsx") %>% 
                   write.xlsx(. , paste0( dir, "side_by_side_comparison/sCD14_genus_20_1/",
                                          "0Sum_sCD14_genus_20_1_sum.xlsx" ) )
#### HIV ####
name = "HIV_genus_20_1"
dir_tmp = paste0("side_by_side_comparison/", name, "/")
nodes_file = "10.1nodes.xlsx"

strong_nodes_diff( nodes = paste0( dir_tmp, nodes_file), 
                   diff_file =  "HIV_Untreated_DE_results.xlsx",
                   corr_file = "HIV_against_features.xlsx" ) %>% 
              write.xlsx(. , paste0( dir, dir_tmp, "0Sum_", name, "_sum.xlsx" ) )

nodes_file = "30.1nodes.xlsx"
strong_nodes_diff( nodes = paste0( dir_tmp, nodes_file), 
                   diff_file =  "HIV_Untreated_DE_results.xlsx",
                   corr_file = "HIV_against_features.xlsx" ) %>% 
  write.xlsx(. , paste0( dir, dir_tmp, "03Sum_", name, "_sum.xlsx" ) )

#### two Omics ####
name = "two_Omics_genus_20_1"
dir_tmp = paste0("side_by_side_comparison/", name, "/")
nodes_file = "40.3nodes.xlsx"

strong_nodes_diff( nodes = paste0( dir_tmp, nodes_file), 
                   diff_file =  "HIV_Untreated_DE_results.xlsx") %>% 
  write.xlsx(. , paste0( dir, dir_tmp, "0Sum_", name, "_sum.xlsx" ) )

#### LTA ####
"Genus_20_1_LTA01nodes_against_Pheno.xlsx"
name = "LTA_genus_20_1"
dir_tmp = paste0("side_by_side_comparison/", name, "/")

nodes_file = "10.2nodes.xlsx"

strong_nodes_diff( nodes = paste0( dir_tmp, nodes_file), 
                   diff_file =  "HIV_Untreated_DE_results.xlsx",
                   corr_file = "LTA_against_features.xlsx" ) %>% 
  write.xlsx(. , paste0( dir, dir_tmp, "0Sum_", name, "_sum.xlsx" ) )

### LPS ##
name = "LPS_genus_20_1"
dir_tmp = paste0("side_by_side_comparison/", name, "/")

nodes_file = "10.2nodes.xlsx"

strong_nodes_diff( nodes = paste0( dir_tmp, nodes_file), 
                   diff_file =  "HIV_Untreated_DE_results.xlsx",
                   corr_file = "LPS_against_features.xlsx"  ) %>% 
  write.xlsx(. , paste0( dir, dir_tmp, "0Sum_", name, "_sum.xlsx" ) )

## sCD14 family ##
name = "sCD14_family_20_1"
dir_tmp = paste0("side_by_side_comparison/", name, "/")

nodes_file = "10.2nodes.xlsx"

strong_nodes_diff( nodes = paste0( dir_tmp, nodes_file), 
                   diff_file =  "HIV_Untreated_DE_results.xlsx") %>% 
  write.xlsx(. , paste0( dir, dir_tmp, "0Sum_", name, "_sum.xlsx" ) )

#### sCD14_genus_40_1_034 ##
name = "sCD14_genus_40_1_034"
dir_tmp = paste0("side_by_side_comparison/", name, "/")

nodes_file = "10.1nodes.xlsx"

strong_nodes_diff( nodes = paste0( dir_tmp, nodes_file), 
                   diff_file =  "HIV_Untreated_DE_results.xlsx") %>% 
  write.xlsx(. , paste0( dir, dir_tmp, "0Sum_", name, "_sum.xlsx" ) )

nodes_file = "20.1nodes.xlsx"

strong_nodes_diff( nodes = paste0( dir_tmp, nodes_file), 
                   diff_file =  "HIV_Untreated_DE_results.xlsx") %>% 
  write.xlsx(. , paste0( dir, dir_tmp, "02Sum_", name, "_sum.xlsx" ) )



######## ##########
## Using sCD14_genus_20, TwoOmics_genus_20, 
 ## HIV_genus_20, LPS_genus_20, LTA_genus_20, 
## sCD14_family_20, sCD14_genus_40_1, 
## sCD14_genus_40_2, sCD14_genus_60_2, sCD14_genus_60_3 as id variables
##########

############ diff Enrich ###########
dir = "~/Documents/gitlab/Omics_Integration/DataProcessed/"

##
## SmCCNet::
