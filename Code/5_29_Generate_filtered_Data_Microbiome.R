########## source then import a function to get filtered microbiome data at that level ######
load_filtered_micro_level <- function(level, prevalence, RA, wd){
  #########
  # level  c("phylum", "Order", "family", "Genus", "Species" )
  # prevalence and RA are positive integers
  # RA = "Ubuntu" "ASUS"
  ##########
  library(openxlsx)
  library(tidyverse)
  options(stringsAsFactors = F)
  if(wd == "ASUS"){
    # SET UP THE folder
    folder =  c('C:/Users/hithr/Documents/Stats/gitlab/tidi_MIBI/')
  }else if (wd == "Ubuntu"){
    
    folder =  c("~/Documents/gitlab/tidi_MIBI/")
  }
  # using the tidi_MIBI functions by Charlie
  source( paste(folder, "Code_charlie/tidi_MIBI.R", sep = "") )
  ############## 26 samples ###########
  small_lib = c("MIHIV998")
  ## shared sample size across datasets
  shared_sam = c( "MIHIV124", "MIHIV132", "MIHIV138", "MIHIV154", "MIHIV178", "MIHIV255", 
                  "MIHIV278", "MIHIV286", "MIHIV323", "MIHIV361", "MIHIV391", "MIHIV404", 
                  "MIHIV428", "MIHIV493", "MIHIV582", "MIHIV594", "MIHIV648", "MIHIV683", 
                  "MIHIV708", "MIHIV716", "MIHIV819", "MIHIV825", "MIHIV839", "MIHIV914", 
                  "MIHIV947", "MIHIV972", "MIHIV998")
  final_sam = shared_sam[ shared_sam %nin% small_lib]
  ########### meta data #############
  #### meta data #########
  old.meta = read.table(paste(folder, "Data/SEQ020_Wilson1_metadata_all_subjects_24Aug2018.txt", 
                               sep = ""), header = T ) 
  sample = gsub( "B", "", old.meta$Lib)
  meta = old.meta %>% dplyr::mutate(Lib = sample) %>% dplyr::filter(Lib %in% final_sam) 
  ######## read microbiome data ##########
  read_data <- function(folder, datafile, header){
    read.table(paste(folder, datafile, sep = ""), header = header)
  }
  read_data <- function(folder, datafile, header){
    read.table(paste(folder, datafile, sep = ""), header = header)
  }
  
  ############## filter the data by sample #############
  #### filter the sample by the shared set of samples
  # the naming convention, a B at the end of the name
  filter_micro <- function(data, samples){
    data[1,-1] = gsub("B", "", data[1,-1])
    name = c("OTU_Name", samples)
    return(data[, (data[1,] %in% name)]  )
  }
  
  ###### having column names#######
  filter_micro_t <- function(data, samples){
    colnames(data)[-1] = gsub("B", "", colnames(data)[-1])
    name = c("OTU_Name", samples)
    return(data[,  colnames(data) %in% name]  )
  }
  
  # load in phylum data with headers
  phy_t <- read_data(folder, "Data/phylum_biopsy_cts_24Aug2018.txt", T)  %>% 
    filter_micro_t(., final_sam)
  ord_t <- read_data(folder, "Data/order_biopsy_cts_24Aug2018.txt", T) %>% 
    filter_micro_t(., final_sam)
  fam_t <- read_data(folder, "Data/family_biopsy_cts_24Aug2018.txt", T) %>% 
    filter_micro_t(., final_sam)
  gen_t <- read_data(folder, "Data/alltaxa_biopsy_cts_24Aug2018.txt", T) %>% 
    filter_micro_t(., final_sam)
  spe_t <- read_data(folder, "Data/species_biopsy_cts_24Aug2018.txt", T) %>% 
    filter_micro_t(., final_sam)
  # first row as the column name 
  header.true <- function(df) {
    names(df) <- as.character(unlist(df[1,]))
    df[-1,]
  }
  ############ The shared sample across datasets #############
  phy <- read_data(folder, "Data/phylum_biopsy_cts_24Aug2018.txt", F) %>% filter_micro(., final_sam)
  ord <- read_data(folder, "Data/order_biopsy_cts_24Aug2018.txt", F) %>% filter_micro(., final_sam)
  fam <- read_data(folder, "Data/family_biopsy_cts_24Aug2018.txt", F) %>% filter_micro(., final_sam)
  gen <- read_data(folder, "Data/alltaxa_biopsy_cts_24Aug2018.txt", F) %>% filter_micro(., final_sam)
  spe <- read_data(folder, "Data/species_biopsy_cts_24Aug2018.txt", F) %>% filter_micro(., final_sam)
  
  ############## filter the data by prevalence and RA #############
  # level c("phylum", "order", "family", "genus", "species" )
  if ( level == "phylum"){
  mibi.set <-  tidi_MIBI(phylum = phy,
                         meta = meta,
                         ## Prevalence cutoff 
                         prev_cutoff = prevalence, 
                         ## Relative abundance cutoff  
                         ra_cutoff = RA,   
                         ## unclassified
                         unc = TRUE)  
  } 
  else if (level == "order"){
    mibi.set <-  tidi_MIBI(class = ord,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff  
                           ra_cutoff = RA,   
                           unc = TRUE) 
  }
  else if (level == "family"){
    mibi.set <-  tidi_MIBI(family = fam,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff  
                           ra_cutoff = RA,   
                           unc = TRUE) 
  }
  else if (level == "genus"){
    mibi.set <-  tidi_MIBI(genus = gen,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff 
                           ra_cutoff = RA,   
                           unc = TRUE) 
  }
  else if (level == "species"){
    mibi.set <-  tidi_MIBI(genus = spe,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff 
                           ra_cutoff = RA,   
                           unc = TRUE) 
  }
  mibi.set$ra <- (mibi.set$cts/ mibi.set$Total)*100
  # long format to wide format, RA
  mibi_filter_ra <-  mibi.set  %>% 
    dplyr::select(Lib, Taxa, ra) %>%
    spread(., Taxa, ra) %>% as.data.frame()
  # long format to wide 
  mibi_filter_clr <-  mibi.set  %>% 
    dplyr::select(Lib, Taxa, clr) %>%
    spread(., Taxa, clr) %>% as.data.frame()
  return(list(mibi_filter_ra, mibi_filter_clr))
}