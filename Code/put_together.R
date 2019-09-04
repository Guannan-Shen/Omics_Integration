library(SmCCNet)
## Load in data 

## Cross validation to find the l1 and l2

########## run SmCCNet wrapper #########
run_SmCCNet <- function(X1, X2, Y, l1, l2, s1, s2, weights, n_na){
  print("weights can be NULL or a length 3 vector")
  runsmccnet <- function(X1, X2, Y, l1, l2, s1, s2, weights){
    x = cbind(X1, X2)
    corr = cor(x, method = "pearson")
    p1 = ncol(X1)
    p2 = ncol(X2)
    n = nrow(X1)
    AbarLabel = c(colnames(cbind(X1, X2)))
    set.seed(123)
    W1 = getRobustPseudoWeights(X1, X2, Trait = Y, Lambda1 = l1,
                                Lambda2 = l2, s1 = s1, s2 = s2, NoTrait = FALSE, FilterByTrait = FALSE,
                                SubsamplingNum = 1000, CCcoef = weights, trace = FALSE)
    abar = getAbar(W1, P1 = p1, FeatureLabel = AbarLabel)
    modules = getMultiOmicsModules(abar, P1 = p1, PlotTree = T, CutHeight = 1 - 0.1^10)
    for(idx in 1:length(modules)){
      print(plotMultiOmicsNetwork(abar, corr, modules, ModuleIdx = idx, P1 = p1, FeatureLabel = AbarLabel,
                                  EdgeCut = 0))
    }
    return(modules)
  }
  if(missing(n_na)){
    runsmccnet(X1, X2, Y, l1, l2, s1, s2, weights)
  }
  else {
    runsmccnet(X1[-n_na, ], X2[-n_na, ], Y[-n_na, ], l1, l2, s1, s2, weights)
  }
}