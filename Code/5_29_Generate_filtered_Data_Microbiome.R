library(compositions)
library(stats)

########## source then import a function to get filtered microbiome data at that level ######

######################### clean microbiome data #####################
## with the above function, the microbiome data has been filtered 
## and the centered log ratio transformation has been applied to the Relative Abundance
## need to scale to mean 0 variance 1
rescale_microbiome <- function(data){
  # the df, eg list[[2]] from the above return has a column Lib
  df = data %>% as.data.frame()  
  results = data.frame(apply(df, 2, scale)) 
  rownames(results) = rownames(df)
  return(results)
}

############# using all 27 samples #########################
load_filtered_micro_level_samples <- function(level, prevalence, RA, wd, collapse){
  #########
  # level  c("phylum", "order", "family", "genus", "species" )
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
  # final_sam = shared_sam[ shared_sam %nin% small_lib]
  final_sam = shared_sam
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
  ######### get library size ################
  lib_s = phy_t[1,-1]
  ################# get RA for all taxa, no filtering, might be useful if we do transformation other than clr ######
  get_RA <- function(data){
    ######## pull out information from data #############
    lib_s = data[1,-1]; samples = colnames(data)[-1]
    taxa = data[-1, 1]; data_core = data[-1,-1]
    #### prevalence ####
    non_zero = apply(data_core , 1, function(x){
      sum(x != 0)
    } ) %>% as.data.frame() %>% .[,1]
    
    ####### calculate RA ############
    data_ra = matrix(NA, nrow = nrow(data_core))
    for( i in 1:ncol(data_core) ){
      # do not round RA
      ra = data_core[,i]/ as.numeric(lib_s[i] ) *100
      data_ra = cbind(data_ra,  ra)
    }
    data_ra = data.frame( data_ra[, -1] )
    row.names( data_ra ) = taxa
    colnames( data_ra ) = samples
    #### avg_RA #####
    avg_ra = apply(data_ra, 1, mean) %>% as.data.frame() %>% .[,1] %>% round(., 4)
    
    res = data_ra  %>% dplyr::mutate(Taxas = taxa,
                                     Prevalence = non_zero,
                                     Avg_RA = avg_ra)  %>% 
      plyr::arrange( plyr::desc(Prevalence), plyr::desc(Avg_RA)) %>% 
      dplyr::mutate(Num = 1:nrow(.)) %>% 
      dplyr::select(Num, Taxas, Prevalence, Avg_RA, everything())
    return(res)
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
    data_ra = get_RA(phy_t)
  } 
  else if (level == "order"){
    mibi.set <-  tidi_MIBI(class = ord,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff  
                           ra_cutoff = RA,   
                           unc = TRUE) 
    data_ra = get_RA(ord_t)
  }
  else if (level == "family"){
    mibi.set <-  tidi_MIBI(family = fam,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff  
                           ra_cutoff = RA,   
                           unc = TRUE) 
    data_ra = get_RA(fam_t)
  }
  else if (level == "genus"){
    mibi.set <-  tidi_MIBI(genus = gen,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff 
                           ra_cutoff = RA,   
                           unc = TRUE) 
    data_ra = get_RA(gen_t)
  }
  else if (level == "species"){
    mibi.set <-  tidi_MIBI(genus = spe,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff 
                           ra_cutoff = RA,   
                           unc = TRUE) 
    data_ra = get_RA(spe_t)
  }
  # microbiome relative abundance should sum to 1 not 100
  mibi.set$ra <- (mibi.set$cts/ mibi.set$Total)*100
  # long format to wide format, RA
  # long format to wide 
  mibi_filter_ra <-  mibi.set  %>% 
    dplyr::select(Lib, Taxa, ra) %>%
    spread(., Taxa, ra) %>% as.data.frame() %>% column_to_rownames("Lib")
  Libs <- rownames(mibi_filter_ra)
  ######### filter genus ##########33
  ### collapse unclassified to the other or not ##########
  if(collapse){
    
  if (level == "genus"){
    non_classified <- c("Cyano:4C0d-2", "Bacte:Bacte:S24-7", 
                        "Prote:Betaproteobacteria", "Tener:Molli:RF9"
                        ## 20% 0%
                        ,"Actin:Actinobacteria", "Cyano:Chloroplast",
                        "Firmi:Bacilli", "Firmi:Clost:Family-XIII-Incertae-Sedis",
                        "Prote:Gammaproteobacteria")
    taxa_all <- colnames(mibi_filter_ra)
    taxa_tofilter <- taxa_all [taxa_all  %in% non_classified]
    ############ combine all to Other with rowSums and mutate ######33
    mibi_filter_ra <- mibi_filter_ra  %>% dplyr::mutate(Other = Other + 
                                                   rowSums(mibi_filter_ra[, taxa_tofilter]%>% as.data.frame()) +
                                                   #### family level #######
                                                   rowSums(dplyr::select(mibi_filter_ra, ends_with("aceae")) ) +
                                                   ###### order level ##########
                                                   rowSums(dplyr::select(mibi_filter_ra, ends_with("ales")) ) ) %>% 
      dplyr::select(-taxa_tofilter, -ends_with("aceae"), -ends_with("ales"))
    rownames(mibi_filter_ra) <- Libs
    ###### library(compositions)
    mibi_filter_clr <- clr( mibi_filter_ra)   %>% as.data.frame()
    rownames(mibi_filter_clr) <- Libs
  }
  else if (level == "family"){
    non_classified <- c("Cyano:4C0d-2", "Firmi:Clost:Finegoldia",
                        "Prote:Betaproteobacteria", "Tener:Molli:RF9"
                        ## 20% 0%
                        ,"Actin:Actinobacteria", "Cyano:Chloroplast",
                        "Firmi:Bacilli","Prote:Gammaproteobacteria",
                        "Firmi:Clost:Peptoniphilus")
    taxa_all <- colnames(mibi_filter_ra)
    taxa_tofilter <- taxa_all [taxa_all  %in% non_classified]
    ############ combine all to Other with rowSums and mutate ######33
    mibi_filter_ra <- mibi_filter_ra  %>% dplyr::mutate(Other = Other + 
                                                          rowSums(mibi_filter_ra[, taxa_tofilter] %>% as.data.frame()) +
                                                          #### family level #######
                                                        # rowSums(dplyr::select(mibi_filter_ra, ends_with("aceae")) ) +
                                                          ###### order level ##########
                                                        rowSums(dplyr::select(mibi_filter_ra, ends_with("ales")) ) ) %>% 
      dplyr::select(-taxa_tofilter, -ends_with("ales"))
    rownames(mibi_filter_ra) <- Libs
    ###### library(compositions)
    mibi_filter_clr <- clr( mibi_filter_ra)   %>% as.data.frame()
    rownames(mibi_filter_clr) <- Libs
  }
  } 
  
  else{ ############ not collapse unclassified to other #############
    if (level == "genus"){
      taxa_all <- colnames(mibi_filter_ra)
      rownames(mibi_filter_ra) <- Libs
      ###### library(compositions)
      mibi_filter_clr <- clr( mibi_filter_ra)   %>% as.data.frame()
      rownames(mibi_filter_clr) <- Libs
    }
    else if (level == "family"){
      taxa_all <- colnames(mibi_filter_ra)
      rownames(mibi_filter_ra) <- Libs
      ###### library(compositions)
      mibi_filter_clr <- clr( mibi_filter_ra)   %>% as.data.frame()
      rownames(mibi_filter_clr) <- Libs
    
    }
  }
  ## return RA and clr, RA of all Taxa and library size 
  return(list(filtered_RA = mibi_filter_ra, 
              filtered_clr = mibi_filter_clr, 
              full_RA = data_ra, 
              library_size = lib_s))
}


#################### get the total number of taxa ################
get_n_taxa <- function(level, prevalence, RA, wd){
  data = load_filtered_micro_level_samples(level,  prevalence, RA, wd)
  df = data[[2]] %>% as.data.frame()
  n = ncol(df)
  return(n)
}

# library(tidy.micro)

################# test run ##############
# prev <- 20
# ra <- 1
# # omcis_name <- "30_2_Global_100_50_Genus"
# 
# micro_data <- load_filtered_micro_level_samples("family",
#                                                 prevalence = prev, RA = ra, wd = "Ubuntu")
# micro_clr <- micro_data[[2]] %>% as.data.frame()
# # rescale to mean 0 and variance 1
# mibi <- rescale_microbiome(micro_clr)
# mibi %>% dim()
# colnames(mibi)
# clrInv(mibi) %>% rowSums()
# mibi[1:5, 1:3]
# test[1:5,1:3]
# test %>% dim()
# sum(rownames(mibi) != clin$ID)
# sum(rownames(mibi) != rownames(filtered_rlog) )
# 
# mibi_ra [1:5, 1:3]
# micro_data[[1]][1:5, 1:3]
# 
# mibi_outlier <- grubbs_df(mibi, 2, 10)$fdr > 0.05
# print(sum(mibi_outlier))
# 
# 
# test <- clr( mibi_ra)   %>% as.data.frame() 
# 
# hist(test[15,]%>% as.numeric())
# hist(micro_clr[15,] %>% as.numeric())

 
############ test data set, without the function ##########
 # folder =  c("~/Documents/gitlab/tidi_MIBI/")
 # read_data <- function(folder, datafile, header){
 #   read.table(paste(folder, datafile, sep = ""), header = header)
 # }
 # 
 # filter_micro <- function(data, samples){
 #   data[1,-1] = gsub("B", "", data[1,-1])
 #   name = c("OTU_Name", samples)
 #   return(data[, (data[1,] %in% name)]  )
 # }
 # ## shared sample size across datasets
 # shared_sam = c( "MIHIV124", "MIHIV132", "MIHIV138", "MIHIV154", "MIHIV178", "MIHIV255", 
 #                 "MIHIV278", "MIHIV286", "MIHIV323", "MIHIV361", "MIHIV391", "MIHIV404", 
 #                 "MIHIV428", "MIHIV493", "MIHIV582", "MIHIV594", "MIHIV648", "MIHIV683", 
 #                 "MIHIV708", "MIHIV716", "MIHIV819", "MIHIV825", "MIHIV839", "MIHIV914", 
 #                 "MIHIV947", "MIHIV972", "MIHIV998")
 # # final_sam = shared_sam[ shared_sam %nin% small_lib]
 # final_sam = shared_sam
 # old.meta = read.table(paste(folder, "Data/SEQ020_Wilson1_metadata_all_subjects_24Aug2018.txt", 
 #                             sep = ""), header = T ) 
 # sample = gsub( "B", "", old.meta$Lib)
 # meta = old.meta %>% dplyr::mutate(Lib = sample) %>% dplyr::filter(Lib %in% final_sam) 
 # gen <- read_data(folder, "Data/alltaxa_biopsy_cts_24Aug2018.txt", F) %>% filter_micro(., final_sam)
 # mibi.set <-  tidi_MIBI(genus = gen,
 #                        meta = meta,
 #                        ## Prevalence cutoff 
 #                        prev_cutoff = 20, 
 #                        ## Relative abundance cutoff 
 #                        ra_cutoff = 1,   
 #                        unc = TRUE) 
 # mibi.set$ra <- (mibi.set$cts/ mibi.set$Total)
 # mibi_filter_ra <- mibi.set  %>% 
 #   dplyr::select(Lib, Taxa, ra) %>%
 #   spread(., Taxa, ra) %>% as.data.frame() %>% column_to_rownames("Lib")
 # Libs <- rownames(mibi_filter_ra)
 # 
 # non_classified <- c("Cyano:4C0d-2", "Bacte:Bacte:S24-7", 
 #                     "Prote:Betaproteobacteria")
 # taxa_all <- colnames(mibi_filter_ra)
 # taxa_tofilter <- taxa_all [taxa_all  %in% non_classified]
 # ############ combine all to Other with rowSums and mutate ######33
 # mibi_filter_ra <- mibi_filter_ra  %>% dplyr::mutate(Other = Other + 
 #                                                       rowSums(mibi_filter_ra[, taxa_tofilter]) +
 #                                                       #### family level #######
 #                                                     rowSums(dplyr::select(mibi_filter_ra, ends_with("aceae")) ) +
 #                                                       ###### order level ##########
 #                                                     rowSums(dplyr::select(mibi_filter_ra, ends_with("ales")) ) ) %>% 
 #   dplyr::select(-taxa_tofilter, -ends_with("aceae"), -ends_with("ales"))
 # rownames(mibi_filter_ra) <- Libs
 # ###### library(compositions)
 # mibi_filter_clr <- clr( mibi_filter_ra)   %>% as.data.frame()
 # rownames(mibi_filter_clr) <- Libs
 # 
 # 

 ##############################################
 #  library(compositions) clr, clrInv, subjects by row
 ##############################################
 
 ########### using compositions clr is better than the original clr of tidy.micro #########



################### remove the smallest library ##################
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
  ######### get library size ################
  lib_s = phy_t[1,-1]
  ################# get RA for all taxa, no filtering, might be useful if we do transformation other than clr ######
  get_RA <- function(data){
    ######## pull out information from data #############
    lib_s = data[1,-1]; samples = colnames(data)[-1]
    taxa = data[-1, 1]; data_core = data[-1,-1]
    #### prevalence ####
    non_zero = apply(data_core , 1, function(x){
      sum(x != 0)
    } ) %>% as.data.frame() %>% .[,1]
    
    ####### calculate RA ############
    data_ra = matrix(NA, nrow = nrow(data_core))
    for( i in 1:ncol(data_core) ){
      # do not round RA
      ra = data_core[,i]/ as.numeric(lib_s[i] ) *100
      data_ra = cbind(data_ra,  ra)
    }
    data_ra = data.frame( data_ra[, -1] )
    row.names( data_ra ) = taxa
    colnames( data_ra ) = samples
    #### avg_RA #####
    avg_ra = apply(data_ra, 1, mean) %>% as.data.frame() %>% .[,1] %>% round(., 4)
    
    res = data_ra  %>% dplyr::mutate(Taxas = taxa,
                                     Prevalence = non_zero,
                                     Avg_RA = avg_ra)  %>% 
      plyr::arrange( plyr::desc(Prevalence), plyr::desc(Avg_RA)) %>% 
      dplyr::mutate(Num = 1:nrow(.)) %>% 
      dplyr::select(Num, Taxas, Prevalence, Avg_RA, everything())
    return(res)
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
    data_ra = get_RA(phy_t)
  } 
  else if (level == "order"){
    mibi.set <-  tidi_MIBI(class = ord,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff  
                           ra_cutoff = RA,   
                           unc = TRUE) 
    data_ra = get_RA(ord_t)
  }
  else if (level == "family"){
    mibi.set <-  tidi_MIBI(family = fam,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff  
                           ra_cutoff = RA,   
                           unc = TRUE) 
    data_ra = get_RA(fam_t)
  }
  else if (level == "genus"){
    mibi.set <-  tidi_MIBI(genus = gen,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff 
                           ra_cutoff = RA,   
                           unc = TRUE) 
    data_ra = get_RA(gen_t)
  }
  else if (level == "species"){
    mibi.set <-  tidi_MIBI(genus = spe,
                           meta = meta,
                           ## Prevalence cutoff 
                           prev_cutoff = prevalence, 
                           ## Relative abundance cutoff 
                           ra_cutoff = RA,   
                           unc = TRUE) 
    data_ra = get_RA(spe_t)
  }
  # microbiome relative abundance should sum to 1 not 100
  mibi.set$ra <- (mibi.set$cts/ mibi.set$Total)*100
  # long format to wide format, RA
  mibi_filter_ra <-  mibi.set  %>% 
    dplyr::select(Lib, Taxa, ra) %>%
    spread(., Taxa, ra) %>% as.data.frame() %>% column_to_rownames("Lib")
  # long format to wide 
  mibi_filter_clr <-  mibi.set  %>% 
    dplyr::select(Lib, Taxa, clr) %>%
    spread(., Taxa, clr) %>% as.data.frame() %>% column_to_rownames("Lib")
  return(list(mibi_filter_ra, mibi_filter_clr, data_ra, lib_s)) }


