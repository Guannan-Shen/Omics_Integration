############# diagnostic plots ###########
## histogram and boxplot with stat_summary ##
'%nin%' <- Negate('%in%')
options(stringsAsFactors = F)
options(dplyr.width = Inf)

library(readxl)
library(tidyverse)
library(magrittr)
library(tools)
# Correlations with significance levels
# library(Hmisc)
setwd("~/Documents/gitlab/Omics_Integration/DataRaw/hiv_infected_un/")
getwd()

## histogram
# the histogram is the distribution of within-subject variance 
# the data has rownames as the subject ID
hist_var <- function(data, title){
  print("Not for Transcriptome!")
  var = apply(data, 1, function(x){
    var(as.vector(as.matrix(x) ))
  })
  p = ggplot(mapping =  aes(x=var)) +
    geom_histogram(color="black", fill="white") +
    labs(caption = title) +
    theme_bw() 
  print(p)
}

## get within-subject variance, by row
get_var <- function(data, label){
  print("Group is the group label")
  print("Not for Transcriptome!")
  var = apply(data, 1, function(x){
    var(as.vector(as.matrix(x) ))
  }) %>% as.data.frame() %>% dplyr::mutate(group = label) %>% 
    dplyr::rename( values = ".") 
  return(var)
}

## get pearson correlation 
get_corr <- function(df1, df2, label){
  print("No missing values")
  print("Group is the group label")
  # pearson correlations
  df = stats::cor(df1, df2, method = "pearson") %>% as.data.frame() %>% stack() %>% 
    as.data.frame()  %>% dplyr::mutate(group = label) %>%
    select(values, group)
  return(df)
}


# boxplots
# https://stackoverflow.com/questions/38020772/annotate-boxplot-in-ggplot2
# the data is a dataframe contains a column of value and anothe column of grouping. 
# the value could be variance and correlations 
box_values_group <- function(data, title){
  # 
  print("A column called group and a column called values")
 p = ggplot(data, aes(x = group, y = values, fill = group )) + 
    geom_boxplot(width=0.3) +
   # grey scale plot
    scale_fill_grey(start = 0.2, end = 0.8) +
   # get the quantiles annotation
    stat_summary(geom="text", fun.y=quantile,
                 aes(label=sprintf("%1.1f", ..y..) ), color = "black",
                 position=position_nudge(x=0.33), size = 4) +
    theme_bw() +
    labs(caption = title) 
 
 print(p)
}


## grouped density plots
# the data is a dataframe contains a column of value and anothe column of grouping. 
# the value could be variance and correlations 
density_values_group <- function(data, title){
  print("A column called group and a column called values")
  p = ggplot(data, aes(x = values, color = group, fill = group)) + 
    # alpha controls the transparency 
    geom_density(alpha = 0.5) +
    # gray scale color
    # grey scale plot
    scale_colour_grey() +
    scale_fill_grey() +
    theme_bw() +
    labs(caption = title) 
  print(p)
}

## check skewness 
check_skew <- function(df, title){
  data = df %>% as.data.frame() %>% stack() %>% 
    as.data.frame()  
  p = ggplot(data, aes(x = values)) + 
    # alpha controls the transparency 
    geom_density(alpha = 0.5, fill = "grey50") +
    # gray scale color
    # grey scale plot
    # scale_colour_grey() +
    # scale_fill_grey() +
    theme_bw() +
    labs(caption = title) 
  print(p)
}

## elbow plot of filtering microbiome 
elbow_prev <- function(taxa_level, prev_vector, RA){
  # source to get the load_filtered_micro_level function to get clr of RA
  dir = "~/Documents/gitlab/Omics_Integration/"
  source( paste0(dir, "Code/5_29_Generate_filtered_Data_Microbiome.R") )
  y_ntaxa = c(NULL)
  for (i in prev_vector){
    data = load_filtered_micro_level_samples(taxa_level,  prevalence = i, RA = RA, wd = "Ubuntu")
    n = data[[2]] %>% as.data.frame() %>% ncol()
    y_ntaxa = c(y_ntaxa, n)
  }
  coords = paste("(", prev_vector, ", ", y_ntaxa, ")", sep="")
  p = ggplot(mapping =  aes(x = prev_vector, y = y_ntaxa )) + 
    theme_bw() +
    geom_line() +
    labs(x = paste("Prevalence %","(", "While relative abundance =",RA ,"%)", sep = " "),
         y = paste0("n taxa at ", toTitleCase(taxa_level), " level") ) + 
    geom_label(aes(prev_vector, y_ntaxa, label=coords))
  print(p)
  ggsave(filename = paste0("~/Documents/gitlab/Omics_Integration/Reports/plots/","Prev", 
               taxa_level, RA, ".tiff"),
         dpi = 300, compression = "lzw" )
}


elbow_RA <- function(taxa_level, prev, RA_vector){
  # source to get the load_filtered_micro_level function to get clr of RA
  dir = "~/Documents/gitlab/Omics_Integration/"
  source( paste0(dir, "Code/5_29_Generate_filtered_Data_Microbiome.R") )
  y_ntaxa = c(NULL)
  for (i in RA_vector){
    data = load_filtered_micro_level_samples(taxa_level,  prevalence = prev, RA = i, wd = "Ubuntu")
    n = data[[2]] %>% as.data.frame() %>% ncol()
    y_ntaxa = c(y_ntaxa, n)
  }
  coords = paste("(", RA_vector, ", ", y_ntaxa, ")", sep="")
  p = ggplot(mapping =  aes(x = RA_vector, y = y_ntaxa )) + 
    theme_bw() +
    geom_line() +
    labs(x = paste("Relative Abundance %","(", "While prevalence =",prev ,"%)", sep = " "),
         y = paste0("n taxa at ", toTitleCase(taxa_level), " level") ) + 
    geom_label(aes(RA_vector, y_ntaxa, label=coords))
  print(p)
  ggsave(filename =  paste0( taxa_level, prev, "RA.tiff"), device = NULL, 
         path = "~/Documents/gitlab/Omics_Integration/Reports/plots/",
         dpi = 300, compression = "lzw")
}
# elbow_RA('genus', 0, c(0, 0.5, 1, 2, 3, 4, 5))
# 
# elbow_prev("genus", seq(0, 70, 10), 0)
# elbow_prev("genus", seq(0, 60, 10), 1)
# elbow_prev("genus", seq(0, 70, 10), 2)

contour_prev_RA <- function(taxa_level, prev_vector, RA_vector){
  
}
