############# diagnostic plots ###########
## histogram and boxplot with stat_summary ##
'%nin%' <- Negate('%in%')
options(stringsAsFactors = F)
options(dplyr.width = Inf)

library(reshape2)
library(readxl)
library(tidyverse)
library(magrittr)
library(tools)
library(wesanderson)
library(extrafont)
library(diffEnrich)
library(UpSetR)

# Correlations with significance levels
# library(Hmisc)
setwd("~/Documents/gitlab/Omics_Integration/DataRaw/hiv_infected_un/")
getwd()


########### HOW TO SAVE #############
save_lzw <- function(title){
  ggsave(filename =  paste0( title, ".tiff"), device = NULL,
         path = "~/Documents/gitlab/Omics_Integration/DataProcessed/plots/non_collapse/",
         dpi = 300, compression = "lzw", 
         width = 10, height = 8, units = "in")
  
}

############# large axis font ##########
  #  theme_bw() +
  # theme(legend.position="bottom", legend.box = "horizontal" ) +
# theme(axis.text.x = element_text(size = 16),
#       axis.text.y = element_text(size = 16),
#       axis.title.x = element_text(size = 18),
#       axis.title.y = element_text(size = 18),
#       legend.text = element_text(size=16),
#       legend.title = element_text(size=16),
#       text=element_text(family="Arial")) 

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
  df = stats::cor(df1, df2, use = "pairwise.complete.obs", 
                  method = "pearson") %>% as.data.frame() %>% stack() %>% 
    as.data.frame()  %>% dplyr::mutate(group = label) %>%
    dplyr::select(values, group)
  return(df)
}

get_corr_1 <- function(df1, label){
  print("No missing values")
  print("Group is the group label")
  # pearson correlations
  df = stats::cor(df1, use = "pairwise.complete.obs", 
                  method = "pearson") %>% as.data.frame() %>% stack() %>% 
    as.data.frame()  %>% dplyr::mutate(group = label) %>%
    dplyr::select(values, group)
  return(df)
}
# get_corr_1(filtered_rlog[, filtered_outlier], "Within-Transcriptome")


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
    theme(legend.position="bottom", legend.box = "horizontal") +
    labs(caption = title) 
 
 print(p)
 ggsave(filename =  paste0( title, "box.tiff"), device = NULL, 
        path = "~/Documents/gitlab/Omics_Integration/Reports/plots/",
        dpi = 300, compression = "lzw")
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
    scale_colour_grey(start = 0.2, end = 0.8) +
    scale_fill_grey(start = 0.2, end = 0.8) +
    theme_bw() +
    theme(legend.position="bottom", legend.box = "horizontal") +
    labs(caption = title) 
  print(p)
  ggsave(filename =  paste0( title, "density.tiff"), device = NULL, 
         path = "~/Documents/gitlab/Omics_Integration/Reports/plots/",
         dpi = 300, compression = "lzw")
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
####### scatter plot to check outliers ##############
dot_groupby <- function(data, y, groupby, dotsize, xlab, ylab){
  p = ggplot(data, aes(x=groupby, y=y)) + 
    # dot plot
    geom_dotplot(binaxis='y', stackdir='center', dotsize = dotsize)+
    scale_colour_grey(start = 0.2, end = 0.8) +
    theme_bw() +
    # mean and 2 std
    stat_summary(fun.data="mean_sdl", fun.args = list(mult=2), 
                 geom="crossbar", width=0.3) +
    labs(x = xlab, y = ylab)
  print(p)
}

# test
# dot_groupby(clin, clin$LPS, clin$Group, 0.8, "HIV Status", "LPS")
# dot_groupby(clin, clin$LPS, 1, 0.8, "The whole cohort", "LPS")
# 
# dot_groupby(clin, clin$CD14, clin$Group, 0.8, "HIV Status", "CD14")
# dot_groupby(clin, clin$CD14, 1, 0.8, "The whole cohort", "CD14")

########### grouped dot plot #############
dot_group_violin <- function(data, groupby, y, xlab, ylab){
  print("Using data$group, data$y for groupby and y arguments")
  p = ggplot(data, aes(x=groupby, y=y)) + 
    # dot plot
    # geom_dotplot(binaxis='y', stackdir='center', dotsize = dotsize)+
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.09) +
    scale_colour_grey(start = 0.2, end = 0.8) +
    # scale_y_continuous(name = waiver(), limits = c(-1.02, 1.02), breaks = c(-1, -0.5, 0, 0.5, 1), 
    #                    labels = c(-1, -0.5, 0, 0.5, 1)) +
    theme_bw() +
    # mean and 2 std
    # stat_summary(fun.data="mean_sdl", fun.args = list(mult=2), 
    #              geom="crossbar", width=0.3) +
    labs(x = xlab, y = ylab) +
    coord_flip() +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17))
  print(p)
}

############ Compare Pearson'r of 3 datasets ###############
compare_corr <- function(mibi, filtered_rlog, Y, level, Y_name){
  micro_name = paste0("Microbiome","(",level, ")")
  # between correlation
  df1 = get_corr(mibi, filtered_rlog, paste0(micro_name,":Transcriptome"))
  df2 = get_corr(Y, mibi, paste0(micro_name,":", Y_name) )
  df3 = get_corr(Y, filtered_rlog, paste0("Transcriptome:", Y_name))
  ## within correlation
    df4 = get_corr_1(mibi, paste0("Within-",micro_name) )
    df5 = get_corr_1(filtered_rlog, "Within-Transcriptome")
    tmp = rbind(df1, df2, df3, df4, df5) %>% as.data.frame()
  ## remove self-correlation
    data = tmp %>% dplyr::filter(values != 1)
    dot_group_violin(data, data$group, data$values, "", "Pearson's r")
}

######## dot plot , boxplot ########
dot_group_box <- function(data, groupby, y, xlab, ylab){
  print("Using data$group, data$y for groupby and y arguments")
  p = ggplot(data, aes(x=groupby, y=y)) + 
    # dot plot
    # geom_dotplot(binaxis='y', stackdir='center', dotsize = dotsize)+
    # geom_violin(trim=FALSE) +
    geom_boxplot(width=0.3) +
    scale_colour_grey(start = 0.2, end = 0.8) +
    # scale_y_continuous(name = waiver(), limits = c(-1.02, 1.02), breaks = c(-1, -0.5, 0, 0.5, 1), 
    #                    labels = c(-1, -0.5, 0, 0.5, 1)) +
    theme_bw() +
    # mean and 2 std
    # stat_summary(fun.data="mean_sdl", fun.args = list(mult=2), 
    #              geom="crossbar", width=0.3) +
    labs(x = xlab, y = ylab) +
    coord_flip() +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17))
  print(p)
}

############# cutoff figures from a df row ###############
micro_cutoffs_prev <- function(prev, taxa_level, ra_list, n_taxa){
  coords = paste(ra_list, "% ", n_taxa, sep="")
  
  p = ggplot(data = df, mapping =  aes(x = ra_list, y = n_taxa )) + 
    theme_bw() +
    geom_line() +
    labs(x = paste0("Relative Abundance Cutoff % ","(", "While prevalence = ",prev ,"%)"),
         y = paste0("N Taxa (", tools::toTitleCase(taxa_level), ")") ) + 
    geom_label(aes(ra_list, n_taxa, label=coords), size = 4.5, label.size = 0.2) +
    scale_x_continuous(limits = c(-0.3, 5.3), breaks = ra_list) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18)) +
    theme(text=element_text(family="Arial"))
  print(p)
  title = paste0(prev, "_prev_Ntaxa_", taxa_level)
  ggsave(filename =  paste0( title, ".tiff"), device = NULL,
         path = "~/Documents/gitlab/Omics_Integration/DataProcessed/plots/non_collapse/",
         dpi = 300, compression = "lzw")
  
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
         y = paste0("n taxa at ", tools::toTitleCase(taxa_level), " level") ) + 
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
         y = paste0("n taxa at ", tools::toTitleCase(taxa_level), " level") ) + 
    geom_label(aes(RA_vector, y_ntaxa, label=coords))
  print(p)
  ggsave(filename =  paste0( taxa_level, prev, "RA.tiff"), device = NULL, 
         path = "~/Documents/gitlab/Omics_Integration/Reports/plots/",
         dpi = 300, compression = "lzw")
}
############### need rerun ############3
# elbow_RA('genus', 0, c(0, 0.5, 1, 2, 3, 4, 5))
# #
# elbow_prev("genus", seq(0, 70, 10), 0)
# 
# # 
# elbow_RA('family', 0, c(0, 0.5, 1, 2, 3, 4, 5))
# #
# elbow_prev("family", seq(0, 70, 10), 0)
# elbow_prev("genus", seq(0, 60, 10), 1)
# elbow_prev("genus", seq(0, 70, 10), 2)
# elbow_RA('genus', 30, c(0, 0.5, 1, 2, 3, 4, 5))
# elbow_RA('genus', 40, c(0, 0.5, 1, 2, 3, 4, 5))


contour_prev_RA <- function(taxa_level, prev_vector, RA_vector){
  
}

## heatmap of features based on correlation 

cor_heatmap <- function(df, method, reorder, hclust_method, text_size, title){
  print("Use pairwise.complete.obs, and methods from pearson, kendall, spearman")
  cormat = cor(df, use = "pairwise.complete.obs", method = method)
  ###### Reordered correlation data visualization #####
  # using correlation as dist
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    # default settings of 
    hc <- hclust(dd, method = hclust_method)
    cormat <-cormat[hc$order, hc$order]
  }
  if(reorder){
    cormat <- reorder_cormat(cormat)
  }else{}
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)] = NA
    return(cormat)
  }
  upper_tri = get_upper_tri(cormat)
  
  # reshape data 
  melted_cormat = melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white", size = 0.25) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name= "Pearson's r") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 16, hjust = 1),
          axis.text.y = element_text(size = 16))+
    coord_fixed() +
    labs(x = "", y = "") + 
    geom_text(aes(Var2, Var1, label = round(value,2) ), color = "black", size = text_size)  +
    theme(
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.5, 0.8),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 10, barheight = 1,
                                 title.position = "top", title.hjust = 0.5)) +
    ######## legend / color bar title ###########
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.text = element_text(size=15),
                legend.title = element_text(size=15),
                text=element_text(family="Arial")  )
  # Print the heatmap
  print(ggheatmap)
  ggsave(filename =  paste0( title, ".tiff"), device = NULL,
         path = "~/Documents/gitlab/Omics_Integration/DataProcessed/plots/non_collapse/",
         dpi = 300, compression = "lzw")
  return(cormat)
}

###### 9_12_Global_rna_integration.Rmd ########33
# clin_pearson <- cor_heatmap(df, "pearson", TRUE, "complete", text_size = 3, title)
# clin_pearson <- cor_heatmap(res, "pearson", TRUE, "complete", text_size = 3, title)
# clin

density_values_ind <- function(data, title){
  print("A column called group and a column called values")
  p = ggplot(data, aes(x = values, color = ind, fill = ind)) + 
    # alpha controls the transparency 
    geom_density(alpha = 0.6) +
    # gray scale color
    # grey scale plot
    scale_colour_grey() +
    # start = 0.2, end = 0.8
    scale_fill_grey() +
    theme_bw() +
    theme(legend.position="bottom", legend.box = "horizontal") +
    labs(caption = title) 
  print(p)
}


###########33 function have better contour plot  ##################
optimize_contour <- function(CVDir, l1_max, l2_max){
  ###### default folder ##########
  top = "~/Documents/gitlab/Omics_Integration/DataProcessed/"
  df = read.csv(paste0(top, CVDir, "ranked_TotalPredictionError.csv"))
  ####### subset ######
  data = df %>% dplyr::select(-X)  %>% dplyr::filter( (l1 <= l1_max) & (l2 <= l2_max)  )
  hmelt = melt(data[ , -3], id.vars = c("l1", "l2")) %>% set_colnames(c("l1", "l2", "Variable",
                                                                        "Error") )
  #####33 format ##########3
  f1 = list(
    family = "Arial, sans-serif",
    size = 31,
    color = "black"
  )
  f2 = list(
    family = "Arial, sans-serif",
    size = 31,
    color = "black"
  )
  
  f3 = list(
    family = "Arial, sans-serif",
    size = 27,
    color = "black"
  )
  f4 = list(
    family = "Arial, sans-serif",
    size = 27,
    color = "black"
  )
  
  a = list(
    title = paste("L1 (Transcriptome)", "(s1 = 0.7",  "s2 = 0.9", ")",sep = " "),
    titlefont = f1,
    showticklabels = TRUE,
    tickfont = f2
  )
  b = list(
    title = "L2 (Microbiome)",
    titlefont = f1,
    showticklabels = TRUE,
    tickfont = f2
  )
  
  c = list(
    title = "Prediction\nError",
    titlefont = f3,
    showticklabels = TRUE,
    tickfont = f4
  )
  ######### contour plot by plot_ly ############3
  contourPlot = plot_ly(hmelt, x = ~l1, y = ~l2, z = ~ Error, type = "contour",
                        colorbar = c ) %>%
    layout(xaxis = a, yaxis = b, showlegend = TRUE)  
  
  # orca not works for me
  export(contourPlot, paste0( CVDir, "refined_L1L2.png"))
  
}

################# test run upset plot ################
######## upset plot ############
df = read.xlsx( paste0("~/Documents/gitlab/Omics_Integration/DataProcessed/side_by_side_comparison/",
                  "all_0.1_nodes_genes.xlsx"))
# example of list input (list of named vectors)
listInput <- list( `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
                   `Two Omics (Genus 20%)` = df$TwoOmics_genus_20 %>% na.omit(),
                   `HIV Status (Genus 20%)` = df$HIV_genus_20 %>% na.omit(),
                   `LPS (Genus 20%)` = df$LPS_genus_20 %>% na.omit(),
                   `LTA (Genus 20%)` = df$LTA_genus_20 %>% na.omit())
## tiff(paste0("~/Documents/gitlab/Omics_Integration/DataProcessed/side_by_side_comparison/",
##            "genes_0.1_nodes_overlapping.tiff"), res = 300)
upset(fromList(listInput), 
      order.by = c("freq", "degree" ), 
      #  keep.order = TRUE,
          point.size = 2,
          line.size = 2,
          text.scale = 2) 


# example of list input (list of named vectors)
listInput <- list( `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
                   `sCD14 (Family 20%)` = df$sCD14_family_20 %>% na.omit())
## tiff(paste0("~/Documents/gitlab/Omics_Integration/DataProcessed/side_by_side_comparison/",
##            "genes_0.1_nodes_overlapping.tiff"), res = 300)
upset(fromList(listInput), 
      order.by = c("freq", "degree" ), 
      #  keep.order = TRUE,
      point.size = 2,
      line.size = 2,
      text.scale = 2) 

# example of list input (list of named vectors)
listInput <- list( `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
                   `sCD14 (Genus 40% Module 1)` = df$sCD14_genus_40_1 %>% na.omit(), 
                   `sCD14 (Genus 40% Module 2)` = df$sCD14_genus_40_2 %>% na.omit(),
                   `sCD14 (Genus 60% Module 2)` = df$sCD14_genus_60_2 %>% na.omit(), 
                   `sCD14 (Genus 60% Module 3)` = df$sCD14_genus_60_3 %>% na.omit()) 
## tiff(paste0("~/Documents/gitlab/Omics_Integration/DataProcessed/side_by_side_comparison/",
##            "genes_0.1_nodes_overlapping.tiff"), res = 300)
upset(fromList(listInput), 
      order.by = c("freq", "degree" ), 
      #  keep.order = TRUE,
      point.size = 2,
      line.size = 2,
      text.scale = 2) 



## dev.off()
## ggsave(filename = paste0( "genes_0.1_nodes_overlapping.tiff"),device = NULL,
##       path = "~/Documents/gitlab/Omics_Integration/DataProcessed/side_by_side_comparison/" , dpi = 300, compression = "lzw", 
##       width = 10, height = 8, units = "in" )

df = read.xlsx( paste0("~/Documents/gitlab/Omics_Integration/DataProcessed/side_by_side_comparison/",
                       "all_0.1_nodes_taxa.xlsx"))
# example of list input (list of named vectors)
listInput <- list( `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
                   `Two Omics (Genus 20%)` = df$TwoOmics_genus_20 %>% na.omit(),
                   `HIV Status (Genus 20%)` = df$HIV_genus_20 %>% na.omit(),
                   `LPS (Genus 20%)` = df$LPS_genus_20 %>% na.omit(),
                   `LTA (Genus 20%)` = df$LTA_genus_20 %>% na.omit() )
upset(fromList(listInput), 
      order.by = c("freq", "degree" ), 
     #  keep.order = TRUE,
          point.size = 2,
          line.size = 2,
          text.scale = 2) 


# example of list input (list of named vectors)
listInput <- list( `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
                   `sCD14 (Family 20%)` = df$sCD14_family_20 %>% na.omit())
## tiff(paste0("~/Documents/gitlab/Omics_Integration/DataProcessed/side_by_side_comparison/",
##            "genes_0.1_nodes_overlapping.tiff"), res = 300)
upset(fromList(listInput), 
      order.by = c("freq", "degree" ), 
      #  keep.order = TRUE,
      point.size = 2,
      line.size = 2,
      text.scale = 2) 


# example of list input (list of named vectors)
listInput <- list( `sCD14 (Genus 20%)` = df$sCD14_genus_20 %>% na.omit(), 
                   `sCD14 (Genus 40% Module 1)` = df$sCD14_genus_40_1 %>% na.omit(), 
                   `sCD14 (Genus 40% Module 2)` = df$sCD14_genus_40_2 %>% na.omit(),
                   `sCD14 (Genus 60% Module 2)` = df$sCD14_genus_60_2 %>% na.omit(), 
                   `sCD14 (Genus 60% Module 3)` = df$sCD14_genus_60_3 %>% na.omit()) 
## tiff(paste0("~/Documents/gitlab/Omics_Integration/DataProcessed/side_by_side_comparison/",
##            "genes_0.1_nodes_overlapping.tiff"), res = 300)
upset(fromList(listInput), 
      order.by = c("freq", "degree" ), 
      #  keep.order = TRUE,
      point.size = 2,
      line.size = 2,
      text.scale = 2) 

