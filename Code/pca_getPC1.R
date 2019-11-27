########## get the first PC ########3

get_pc1_listfeature <- function(features_list){
  n = length(features_list)
  pc1_s = vector("list", n)
  for (i in 1:n){
    n_col = ncol(features_list[[i]])
    if (n_col != 0){
      data = features_list[[i]]
      # observation id 
      pid = rownames(data)
      # variable names
      var_names = colnames(data)
      pca = prcomp(data, center = TRUE, scale. = TRUE, retx = TRUE)
      pc1_s[[i]] = pca$x[ ,"PC1"]
    }
  }
  return(pc1_s)
}
## get_pc1_listfeature(features_cut[[1]])
#### then cor.test ########
cor_test_df <- function(df, y){
  mat = apply(df, 2, function(x){
    test = cor.test( x , y,  alternative = "two.sided", method = "pearson")
    res = c(pvalue = test$p.value, corr = test$estimate)
    return(res)
  })
  return(mat)
}

## wilcox ##
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

######### for HIV ########
## wilcox_test_df (features_01_tmp, HIV[,1])

# features_01_tmp = features_cut [[1]][[3]] %>% as.data.frame()
# tmp_cor = wilcox_test_df(features_01_tmp, Y[,1]) %>% round(., 4) %>% as.data.frame() %>% 
#                 t() %>% as.data.frame() 
# id = sim_micro_names(rownames(tmp_cor)) 
# final01_cor <-  tmp_cor  %>% dplyr::mutate(`FDR (correlation)` = 
#                                              round( p.adjust(pvalue, method = "BH"), 4),
#                                            id = id) %>%
#   plyr::arrange(pvalue) %>% dplyr::rename(`p (correlation)` = pvalue) %>%
#   dplyr::select(id, `FDR (correlation)`, everything() )
# rownames(final01_cor) <- id
# write.xlsx( final01_cor, 
#             file = paste0(dir, run, "301nodes_against_Pheno.xlsx"))

######### PLOT FOR PC1 against HIV ###########

library(extrafont)

df_bar = data.frame(gene = mat[,1], 
                    HIV = ifelse(HIV[,1] == 1, "HIV-1-infected", "Uninfected") )
ggplot(data = df_bar, aes(x = HIV, y = gene, fill = HIV)) +
  geom_boxplot(width=0.3) +
  # grey scale plot
  scale_fill_grey(start = 0.4, end = 0.8) +
  theme_bw() +
  theme(legend.position="bottom", legend.box = "horizontal") +
  labs(fill = "HIV Status", x = "HIV Status", y = "PC1") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        text=element_text(family="Arial"))

ggsave(filename = paste0(run, "1PC1_bar.tiff"),device = NULL,
       path = dir, dpi = 300, compression = "lzw" )

df_bar = data.frame(gene = mat[,2], 
                    HIV = ifelse(HIV[,1] == 1, "HIV-1-infected", "Uninfected") )
ggplot(data = df_bar, aes(x = HIV, y = gene, fill = HIV)) +
  geom_boxplot(width=0.3) +
  # grey scale plot
  scale_fill_grey(start = 0.4, end = 0.8) +
  theme_bw() +
  theme(legend.position="bottom", legend.box = "horizontal") +
  labs(fill = "HIV Status", x = "HIV Status", y = "PC1") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        text=element_text(family="Arial"))
ggsave(filename = paste0(run, "3PC1_bar.tiff"),device = NULL,
       path = dir, dpi = 300, compression = "lzw" )