# directory, Ubuntu 
dir = "~/Documents/gitlab/Omics_Integration/"
getwd()
source( paste0(dir, "Code/CV_caret_smcca.R") )

library(SmCCNet)
## Load in data 

## Cross validation to find the l1 and l2

########## run SmCCNet wrapper #########
run_SmCCNet <- function(X1, X2, Y, l1, l2, s1, s2, weights, n_na, NoTrait, EdgeCut){
  print("weights can be NULL or a length 3 vector")
  runsmccnet <- function(X1, X2, Y, l1, l2, s1, s2, weights, NoTrait, EdgeCut){
    x = cbind(X1, X2)
    corr = cor(x, method = "pearson")
    p1 = ncol(X1)
    p2 = ncol(X2)
    n = nrow(X1)
    AbarLabel = c(colnames(cbind(X1, X2)))
    set.seed(123)
    Ws = getRobustPseudoWeights(X1, X2, Y, Lambda1 = l1,
                                Lambda2 = l2, s1 = s1, s2 = s2, NoTrait = NoTrait, FilterByTrait = FALSE,
                                SubsamplingNum = 1000, CCcoef = weights, trace = FALSE)
   
    abar = getAbar(Ws, P1 = p1, FeatureLabel = AbarLabel)
    
    modules = getMultiOmicsModules(abar, P1 = p1, PlotTree = T, CutHeight = 1 - 0.1^10)
    save(Ws, abar, modules, file = paste0(
                                          # CVDir, 
                                          "SmCCNetWeights.RData"))
    for(idx in 1:length(modules)){
      print(plotMultiOmicsNetwork(abar, corr, modules, ModuleIdx = idx, P1 = p1, FeatureLabel = AbarLabel,
                                  EdgeCut = EdgeCut))
    }
    return(modules)
  }
  if(missing(n_na)){
    runsmccnet(X1, X2, Y, l1, l2, s1, s2, weights, NoTrait, EdgeCut)
  }
  else {
    runsmccnet(X1[-n_na, ], X2[-n_na, ], Y[-n_na, ], l1, l2, s1, s2, weights, NoTrait, EdgeCut)
  }
}