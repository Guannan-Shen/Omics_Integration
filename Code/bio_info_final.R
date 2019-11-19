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

############# read in nodes of the most concise networks ##########
strong_nodes_diff <- function(nodes, diff_file){
  
  dir = "~/Documents/gitlab/Omics_Integration/DataProcessed/"
  diff_file = read.xlsx( paste0(dir, "HIV_Untreated_DE_results.xlsx") )
}

dir = "~/Documents/gitlab/Omics_Integration/DataProcessed/"
diff_file = read.xlsx( paste0(dir, "HIV_Untreated_DE_results.xlsx") )
diff_file %>% plyr::arrange(FDR) %>% 
          dplyr::mutate(`FC (log)` = round(logFC, 2),
                        `p adjusted (FDR)` = format.pval(FDR, digits = 3) ) %>% 
          dplyr::select(Symbol, `FC (log)`, `p adjusted (FDR)`) %>% head()

df = read.xlsx( paste0(dir, "side_by_side_comparison/sCD14_genus_20_1/50.4nodes.xlsx") )
