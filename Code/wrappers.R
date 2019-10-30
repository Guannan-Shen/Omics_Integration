# Wrapper to generate Density plot of 3 groups 
wrapper_density <- function(data1, data2, data3, name1, name2, name3, n_na){
  # heads up
  base::print("This is to empirically define the weights for weighted SmCCNet.")
  # missing data position
  print( 
    ifelse( sum(n_na)== 0, 'Index of missing data not specified', paste('n_na =', n_na) ) 
         )
  # generate correlation matrices
  # diagnostic plots and tables
  # directory, Ubuntu 
  dir = "~/Documents/gitlab/Omics_Integration/"
  source( paste0(dir, "Code/ref_plots.R") )
  makeplots <- function(data1, data2, data3, name1, name2, name3){
    df1 = get_corr(data1, data2, paste0(name1, ":", name2))
    # stats::cor(clin$LPS, mibi, use = "pairwise.complete.obs", method = "pearson")
    df2 = get_corr(data1, data3, paste0(name1, ":", name3))
    df3 = get_corr(data2, data3, paste0(name2, ":", name3))
    # plots
    data = rbind(df1, df2, df3)
    box_values_group(data, 
                     paste("Pearson Correlations Summary:", name1, name2, name3) )
    density_values_group(data, 
                         paste("Pearson Correlations Summary:", name1, name2, name3))
  }
  if(!(sum(n_na)== 0)){
    makeplots(data1[-n_na, ], data2[-n_na, ], data3[-n_na, ], name1, name2, name3)
  }
  else{
    makeplots(data1, data2, data3, name1, name2, name3)
  }
}

######## test #######
# wrapper_density(isgs_rlog, mibi, LPS, 
#                 "Core-ISGs", "Genus (Microbiome)", "LPS", n_na = which(is.na( LPS)))
# wrapper_density(isgs_rlog, mibi, CD14, 
#                 "Core-ISGs", "Genus (Microbiome)", "CD14", n_na = which(is.na( CD14)))


wrapper_load_Omics <- function(genelist, micro_level, prev_cutoff, ra_cutoff){
  options(stringsAsFactors = F)
  options(dplyr.width = Inf)
  dir = "~/Documents/gitlab/Omics_Integration/"
  if(genelist == "Inteferome"){
    
  }
  else if(genelist == "Global"){
    
  }
}






