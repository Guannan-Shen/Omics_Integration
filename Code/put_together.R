# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
getwd()
source( paste0(dir, "Code/CV_caret_smcca.R") )

library(SmCCNet)
## Load in data 

## Cross validation to find the l1 and l2

########## run SmCCNet wrapper #########
# This approach can be useful when the phenotype is not quantitative. 
# To choose SsCCA, set NoTrait =  FALSE and FilterByTrait = TRUE
run_SmCCNet <- function(X1, X2, Y, l1, l2, s1, s2, weights, n_na, NoTrait, EdgeCut, bytrait){
  print("weights can be NULL or a length 3 vector")
  runsmccnet <- function(X1, X2, Y, l1, l2, s1, s2, weights, NoTrait, EdgeCut, bytrait){
    x = cbind(X1, X2)
    corr = cor(x, method = "pearson")
    p1 = ncol(X1)
    p2 = ncol(X2)
    n = nrow(X1)
    AbarLabel = c(colnames(cbind(X1, X2)))
    set.seed(123)
    Ws = getRobustPseudoWeights(X1, X2, Y, Lambda1 = l1,
                                Lambda2 = l2, s1 = s1, s2 = s2, NoTrait = NoTrait, FilterByTrait = bytrait,
                                SubsamplingNum = 1000, CCcoef = weights, trace = FALSE)
   
    abar = getAbar(Ws, P1 = p1, FeatureLabel = AbarLabel)
    
    modules = getMultiOmicsModules(abar, P1 = p1, PlotTree = T, CutHeight = 1 - 0.1^10)
    save(Ws, abar, modules, file = paste0(
                                          CVDir, 
                                          "SmCCNetWeights.RData"))
    for(idx in 1:length(modules)){
      print(plotMultiOmicsNetwork(abar, corr, modules, ModuleIdx = idx, P1 = p1, FeatureLabel = AbarLabel,
                                  EdgeCut = EdgeCut))
    }
    return(modules)
  }
  if(missing(n_na)){
    runsmccnet(X1, X2, Y, l1, l2, s1, s2, weights, NoTrait, EdgeCut, bytrait)
  }
  else {
    runsmccnet(X1[-n_na, ], X2[-n_na, ], Y[-n_na, ], l1, l2, s1, s2, weights, NoTrait, EdgeCut, bytrait)
  }
}


edgecut_by <- function(X1, X2, edgeCut){
  print("Load proper similarity matrix (abar) and modules!")
  x = cbind(X1, X2)
  corr = cor(x, method = "pearson")
  p1 = ncol(X1)
  p2 = ncol(X2)
  n = nrow(X1)
  AbarLabel = c(colnames(cbind(X1, X2)))
  EdgeCut = edgeCut
  
  for(idx in 1:length(modules)){
    print(
      plotMultiOmicsNetwork(abar, corr, modules, ModuleIdx = idx, P1 = p1, FeatureLabel = AbarLabel,
                            EdgeCut = EdgeCut))
  }
}

sim_micro_names <- function(micro_names){
  p2 = length(micro_names)
  micro_names_sim = rep(NULL, p2)
  
  for (i in 1:p2){
    new_vector = unlist(strsplit(micro_names[i], fixed = T, split = "."))
    test_name = new_vector[-c(1,2)] %>% paste(., collapse = ".")
    if (length(new_vector) == 1 ) {
      micro_names_sim[i] = micro_names[i]
    } else if ( length(new_vector) == 2 ){  
      micro_names_sim[i] <-  new_vector[-1] %>%
        paste(., collapse = ".")
    } else if(test_name == 2){
      micro_names_sim[i] <- "4C0d-2"
    }
      else{
      micro_names_sim[i] <-  new_vector[-c(1,2)] %>%
        paste(., collapse = ".")}
  }
  return(micro_names_sim)
}



# load(paste0(dir, "SmCCNetWeights.RData"))
# edgecut_by(filtered_rlog[, filtered_outlier], mibi[, mibi_outlier], 0.1)