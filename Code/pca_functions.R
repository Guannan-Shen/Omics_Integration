#############################################################
######### PCA analysis by prcomp is preferred #######
### prcomp is preferred for numerical accuracy #########
###########################################################
## also 'princomp' can only be used with more units than variables


## load microbiome 
########## source then import a function to get filtered microbiome data at that level ######
load_filtered_micro_level_samples <- function(level, prevalence, RA, wd){
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
      arrange( desc(Prevalence), desc(Avg_RA)) %>% 
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
  mibi.set$ra <- (mibi.set$cts/ mibi.set$Total)*100
  # long format to wide format, RA
  mibi_filter_ra <-  mibi.set  %>% 
    dplyr::select(Lib, Taxa, ra) %>%
    spread(., Taxa, ra) %>% as.data.frame() %>% column_to_rownames("Lib")
  # long format to wide 
  mibi_filter_clr <-  mibi.set  %>% 
    dplyr::select(Lib, Taxa, clr) %>%
    spread(., Taxa, clr) %>% as.data.frame() %>% column_to_rownames("Lib")
  return(list(mibi_filter_ra, mibi_filter_clr, data_ra, lib_s))
}

############## PCA analysis of microbiome clr, using plots #######################
pca_microbiome_plots <- function(levels, prevalence, RA, wd){
  #
  print("Since at Phylum level, n_taxa < n_Lib. Please not use Phylum level data.")
  # get the clr data for unfiltered data
  data = load_filtered_micro_level_samples(levels,  prevalence, RA, wd)[[2]]
  # variance plots and 2d pc plot
  pca_all = prcomp(data, center = TRUE, scale. = TRUE, retx = TRUE)
  ## from here, calculate the percentage of variance holded by each pcs
  pca_all_sd = ( pca_all$sdev/sum(pca_all$sdev) )* 100
  ########## pc names ##############
  pc = colnames(pca_all$x)
  # Horizontal bar plot
  pc_sd_df = data.frame(PCs = pc,
                        Variance = pca_all_sd)
  
  p_var = ggplot(data=pc_sd_df , aes(x=1:length(pc), y=Variance)) +
    geom_bar(stat="identity", fill = "Black") +
    ylab("% Variance") +
    scale_x_continuous(name = "Principal components", breaks = 1:length(pc),
                       labels = pc)  +
    # Horizontal bar plot
    coord_flip() +
    theme_minimal()
  print(p_var)
  
  
  print("Check MIHIV998")
  pca12 = ggplot(data = data.frame(pca_all$x), aes(x = PC1, y = PC2, label = rownames(pca_all$x))) +
    geom_point(shape = 19) +
    xlab(paste("PC1 (", round(pc_sd_df$Variance[1],2) ,"%", ")", 
               sep = "") ) +
    ylab(paste("PC2 (", round(pc_sd_df$Variance[2],2) ,"%", ")", 
               sep = "") ) +
    labs(caption = levels) +
    stat_ellipse(type = "t", linetype = "solid") +
    geom_label() +
    theme_bw()
  print(pca12)
  
  pca13 = ggplot(data = data.frame(pca_all$x), aes(x = PC1, y = PC3, label = rownames(pca_all$x))) +
    geom_point(shape = 19) +
    xlab(paste("PC1 (", round(pc_sd_df$Variance[1],2) ,"%", ")", 
               sep = "") ) +
    ylab(paste("PC3 (", round(pc_sd_df$Variance[3],2) ,"%", ")", 
               sep = "") ) +
    labs(caption = levels) +
    stat_ellipse(type = "t", linetype = "solid") +
    geom_label() +
    theme_bw()
  print(pca13)
  
  pca23 = ggplot(data = data.frame(pca_all$x), aes(x = PC2, y = PC3, label = rownames(pca_all$x))) +
    geom_point(shape = 19) +
    xlab(paste("PC2 (", round(pc_sd_df$Variance[2],2) ,"%", ")", 
               sep = "") ) +
    ylab(paste("PC3 (", round(pc_sd_df$Variance[3],2) ,"%", ")", 
               sep = "") ) +
    labs(caption = levels) +
    stat_ellipse(type = "t", linetype = "solid") +
    geom_label() +
    theme_bw()
  print(pca23)
  
}





