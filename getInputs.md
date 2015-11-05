Before get started, these instructions are used to save the data files that are iniatially in .csv format. 

```R
  setwd("/mnt/research/quantgen/ARVALIS2/PIPELINE2015")
  X <- as.matrix(read.csv2("data/X.csv",header=T,row.names=1))
  Y <- read.csv2("data/Y.csv",header=T)
  W <- as.matrix(read.csv2("data/W.csv",header=F,skip=1))
  colnames(W) <- scan("data/W.csv",encoding="latin1",nlines=1,sep=";",what="character")

  mode(W) <- "numeric"
  mode(X) <- "numeric"
  X <- t(X)

  save(X,file="X.RData")
  save(W,file="W.RData")
  save(Y,file="Y.RData")
```

This code is for preparing data: matching genotypic and genotypic information, applying a quality control to genotypes.

```R
## Parameters
  genoFile='/mnt/research/quantgen/ARVALIS2/PIPELINE2015/data/X.RData'
  envCovFile='/mnt/research/quantgen/ARVALIS2/PIPELINE2015/data/W.RData' 
  phenoFile='/mnt/research/quantgen/ARVALIS2/PIPELINE2015/data/Y.RData' 
  outputFolder='/mnt/research/quantgen/ARVALIS2/PIPELINE2015/input/'
  
  thresholdMAF <- 0.05
  thresholdNAfreq <- 0.3
  valid_markers  <-  c(0,1,2)
##

  dir.create(outputFolder) 
  setwd(outputFolder)
  load(phenoFile)
  load(genoFile)
  load(envCovFile)

 ### Matching Y, X, and W ###################################
   IDs=intersect(Y$VAR,rownames(X))
   index <- Y$VAR%in%IDs
   Y=Y[index,]
   W=W[index,]
   X=X[IDs,];   stopifnot(all(Y$VAR%in%rownames(X))); stopifnot(all(rownames(X)%in%Y$VAR)) 
   Y$y <- Y$rdt

 ### Applying quality control to Genotypes ###################################
   freqNA <- apply(X,2,function(x) sum(!x%in%valid_markers)/length(x)) 
   MAF <- apply(X,2,function(x) mean(x[x%in%valid_markers])/2)
   MAF <- ifelse(MAF>0.5,1-MAF,MAF)

  # Remove markers with MAF and freqNA smaller than thresholds
   index <- which(MAF<thresholdMAF | freqNA>thresholdNAfreq)
   if(length(index)>0) X <- X[,-index]

  # Genomic relationship matrix
   Z <- scale(X,center=TRUE,scale=TRUE)
   G <- tcrossprod(Z)/ncol(Z)

 ## Saving data
  save(X,Y,W,G,file="standardized_data.RData")

 ### Creating inputs for models
  # Model 1: Additive model (EL, without markers)
  ETA <- list(ENV=list(~factor(Y$ENV)-1,model='BRR'),
           VAR=list(~factor(Y$VAR)-1,model='BRR')
          )
  save(ETA,file='ETA_EL.RData')
  
  # Model 2: Additive model (EG, with markers)
  myIndex <- match(Y$VAR,rownames(G))
  EVD <- eigen(G)
  PC <- EVD$vectors[,EVD$values>1e-5]
  for(i in 1:ncol(PC)) PC[,i] <- PC[,i]*sqrt(EVD$values[i])
  PC <- PC[myIndex,]
  ETA$VAR <- list(X=PC,model='BRR')
  save(ETA,file='ETA_EG.RData')
  
  # Model 3: Additive model (EGW, with markers and env. cov.)
  W <- scale(W)/sqrt(ncol(W))
  ETA$COV=list(X=W, model='BRR')
  save(ETA,file='ETA_EGW.RData')
  
  # Model 4: GxW Model 1 (EGW+GxW, Model 3 + interactions between markers and env. covariates)
  Omega <- tcrossprod(W)
  Omega <- Omega/mean(diag(Omega))
  GxW <- G[myIndex,myIndex]*Omega
  EVD <- eigen(GxW)
  PC <- EVD$vectors[,EVD$values>1e-5]
  for(i in 1:ncol(PC)) PC[,i] <- PC[,i]*sqrt(EVD$values[i])
  ETA$GxW=list(X=PC, model='BRR')
  save(ETA,file='ETA_EGW_GxW.RData')
  
```
[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
