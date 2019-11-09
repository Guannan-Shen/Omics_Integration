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

# Correlations with significance levels
# library(Hmisc)
setwd("~/Documents/gitlab/Omics_Integration/DataProcessed/")
getwd()


############ tidy the differential HIV untreated datasets #########
df <- read.csv("res.edger.csv") %>% dplyr::select(Gene_ID, Symbol, everything())
write.xlsx(df, "HIV_Untreated_DE_results.xlsx")
