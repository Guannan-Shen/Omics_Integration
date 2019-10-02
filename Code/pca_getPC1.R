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



