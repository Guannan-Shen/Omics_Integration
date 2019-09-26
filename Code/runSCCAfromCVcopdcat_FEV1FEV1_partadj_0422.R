#######################################################################
## Functions

source("/Users/wz/COPD_SCCA/SCCA/Code/SCCAdiagTools.R")
source("/Users/wz/COPD_SCCA/SCCA/Code/HeritSCCASource.R")
source("/Users/wz/COPD_SCCA/SCCA/Code/PMAmodified.R")

### Obtain cross validation (CV) results.


plotCVcontour <- function(CVDir, SCCAmethod, NumSubsamp = 1000){
  plotD <- paste0(CVDir, "Figures/")
  saveD <- paste0(CVDir, "Results/")
  dir.create(plotD); dir.create(saveD)
  
  K <- length(Sys.glob(paste0(CVDir, "*/"))) - 2
  dCor1 <- dCor2 <- NULL
  for(j in 1:K){
    dCorT <- read.csv(paste0(CVDir, "CV_", j, "/", SCCAmethod, 
                             "/SCCA_", NumSubsamp, "_allDeltaCor.csv"))
    dCor1 <- cbind(dCor1, dCorT$CC.Pred..Error)
    #dCor2 <- cbind(dCor2, dCorT$DeltaCor.alt)
  }
  S1 <- rowMeans(dCor1)
  #S2 <- rowMeans(dCor2)
  print(which(S1 == min(S1)))
  #print(which(S2 == min(S2)))
  T12 <- dCorT[ , -1]; T12[ , 3] <- S1; #T12[ , 4] <- S2
  write.csv(T12, file = paste0(saveD, SCCAmethod, "CVmeanDeltaCors.csv"))
  
  PlotEntropyContour(T12[ , 1:3], 
                     paste0(plotD, SCCAmethod, "meanDeltaCorContour.pdf"))
  # PlotEntropyContour(T12[ , c(1:2, 4)],
  #                    paste0(plotD, SCCAmethod, "meanDeltaCorAltContour.pdf"))
}



#######################################################################
# Necessary variables.
# 
# phenoDir <- "copdGeneMicroRna/SCCA/RUVr_mirnak3_mrnak3_GenderRaceAgeSmokCBCFEV1/pctEmph/"
# dataF <- paste0(phenoDir, "Data.Rdata")
#########################################################################

phenoDir <- paste0("~/COPD_SCCA/remove_outliers/FEV1FEV1_partadj_vst_noOut/")
load(paste0(phenoDir, "FEV1FEV1_partadj_vst_noOut.Rdata"))


set.seed(12345)

keep <- which(!is.na(Y.pheno$FEV1pp_utah_P2))
Y <- matrix(Y.pheno$FEV1pp_utah_P2[keep], ncol = 1)

# keep <- which(!is.na(Y.pheno$pctEmph_Thirona_P2))
# Y <- matrix(Y.pheno$pctEmph_Thirona_P2[keep], ncol = 1)

X1 <- t(X.gene[ , keep]); X2 <- t(X.mirna[ , keep])
N <- length(keep); p1 <- ncol(X1); p2 <- ncol(X2)
K <- 4

#dir.create(paste0(resultDir, K, "foldCVpctEmph/" ))

s1 <- 0.7; s2 <- 0.9

pen <- seq(.01, .5, by = .02)
P1P2 <- expand.grid(pen, pen)


CVDir <- paste0(phenoDir, K, "foldCVpctEmph/")
plotD <- paste0(CVDir, "Figures/")
saveD <- paste0(CVDir, "Results/")

# Load data file.
# load(dataF)

bigCor <- cor(cbind(X1, X2, Y))

p1 <- ncol(X1)
AbarLabel <- c(colnames(cbind(X1, X2, Y)))

genenames <- read.csv("/Users/wz/COPD_SCCA/proteinCodingGeneNames.csv")
X1_gene <- as.data.frame(colnames(X1))
X1_gene$id  <- 1:nrow(X1_gene)
X1_gene <- merge(X1_gene, genenames, by.x = "colnames(X1)", by.y = "ensembl_gene_id", all.x = TRUE)
X1_gene <- X1_gene[order(X1_gene$id), ]

AbarLabel <- c(as.character(X1_gene$external_gene_name), colnames(cbind(X2, Y)))


rownames(bigCor) <- colnames(bigCor) <- AbarLabel


# Create contour plots for canonical correlation prediction error. 
subSamp <- 1000

plotCVcontour(CVDir, "SmCCA", NumSubsamp = subSamp)
#plotCVcontour(CVDir, "SsCCA", NumSubsamp = subSamp)


# Run SCCA on the entired dataset based on CV proposed penalty parameters.
for(Method in c("SmCCA")){
  T12 <- read.csv(paste0(CVDir, "Results/", Method, "CVmeanDeltaCors.csv"))[ , -1]
  
  if(Method == "SmCCA"){
    #pen <- which(T12$DeltaCor.alt == min(T12$DeltaCor.alt))
    pen <- which(T12$CC.Pred..Error == min(T12$CC.Pred..Error))
    FilterByTrait <- FALSE
  }else if(Method == "SsCCA"){
    pen <- which(T12$DeltaCor == min(T12$DeltaCor))
    FilterByTrait <- TRUE
  }
  
  Y.scaled <- scale(Y, center = TRUE, scale = TRUE)
  l1 <- T12$l1[pen]; l2 <- T12$l2[pen]
  Ws <- myprocB(X1, X2, Y.scaled, l1, l2, s1, s2, NoTrait = FALSE,
                FilterByTrait = FilterByTrait, Bipartite = TRUE,
                SubsamplingNum = subSamp)
  Ws <- Matrix(Ws)
  Abar <- computeEdgeWeightMatrix(Ws)
  
  save(l1, l2, X1, X2, Y, Y.scaled, s1, s2, Ws, Abar,
       file = paste0(saveD, Method, K, "foldSamp", subSamp, "_", pen, ".Rdata"))
  
  
  sccafile <- paste0(saveD, Method, K, "foldSamp", subSamp, "_", pen, ".Rdata")
  saveresult <- paste0(saveD, Method, K, "foldSamp", subSamp, "_", pen, "Result.Rdata")
  penDir <- paste0(plotD, Method, "pen", pen, "/")
  dir.create(penDir)
  hcPlot <- paste0(penDir, Method, "_", subSamp, "_", pen, "HC.pdf")
  
  R <- GetFuncGrp(sccafile, saveresult, hcPlot, Cutoff = 1,
                  GetPerformance = FALSE)
  
  load(saveresult)
  
  modules <- Leaves2Module(lower.leaves, cutree(hc, h = 1-.1^10))
  mirGeneModule <- lapply(modules, function(s){
    s.min <- min(s)
    s.max <- max(s)
    if(s.min <= p1 & s.max > p1)return(s)
  })
  
  if(length(mirGeneModule) > 1){
    nullSet <- which(sapply(mirGeneModule, is.null))
    if(length(nullSet) > 0){
      mirGeneModule <- mirGeneModule[-nullSet]
    }
  }
  save(pen, penDir, hc, l1, l2, Cutoff, lower.leaves, mirGeneModule,
       Signal, precision, recall,
       file = saveresult)
  
  
  for(edgeCut in c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)){
    netD <- paste0(penDir, "NetworkEdgeCutpt", substring(edgeCut, 3), "/")
    dir.create(netD)
    Abar <- as.matrix(Abar)
    colnames(Abar) <- rownames(Abar) <- AbarLabel[1:nrow(Abar)]
    
    lapply(1:length(mirGeneModule), function(z){
      saveplot.z <- paste0(netD, "SCCA_", subSamp, "_", pen, "_net", z, ".pdf")
      plottitle <- paste0("SCCA_", subSamp, "_", pen, "_net", z)
      PlotReducedNetwork(Abar, bigCor, mirGeneModule, 
                         z, p1 = p1,
                         AddCorrSign = TRUE,
                         SaveFile = saveplot.z, EdgeCut = edgeCut,
                         WeightedEdge = TRUE, ShowNodes = FALSE,
                         PlotTitle = plottitle)
    })
  }
} 

