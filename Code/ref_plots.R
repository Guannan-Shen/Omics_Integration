############# diagnostic plots ###########
## histogram and boxplot with stat_summary ##
'%nin%' <- Negate('%in%')
options(stringsAsFactors = F)

library(readxl)
library(tidyverse)
library(magrittr)
# Correlations with significance levels
# library(Hmisc)
setwd("~/Documents/gitlab/Omics_Integration/DataRaw/hiv_infected_un/")
getwd()

## histogram
# the histogram is the distribution of within-subject variance 
# the data has rownames as the subject ID
hist_var <- function(data, title){
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

