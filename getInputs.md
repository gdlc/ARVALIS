This code is for preparing data: matching genotypic and genotypic information, applying a quality control to genotypes.

```R
## Parameters
  genoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/X_2012_2014.rda'
  envCovFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/W_No_ctr_std.rda' 
  phenoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/Y.rda' 
  outputFolder='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/'
  
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
   W=W_No_ctr_std[index,]
   X=X[IDs,];   stopifnot(all(Y$VAR%in%rownames(X))); stopifnot(all(rownames(X)%in%Y$VAR)) 

 ### Applying quality control to Gentoypes ###################################
   freqNA <- apply(X,2,function(x) sum(!x%in%valid_markers)/length(x)) 
   MAF <- apply(X,2,function(x) mean(x[x%in%valid_markers])/2)
   MAF <- ifelse(MAF>0.5,1-MAF,MAF)

  # Remove markers with MAF and freqNA smaller than thresholds
   index <- which(MAF<thresholdMAF | freqNA>thresholdNAfreq)
   if(length(index)>0) X <- X[,-index]

  # Genomic relationship matrix
   Z <- scale(X,center=TRUE,scale=TRUE)
   G <- tcrossprod(Z)/ncol(Z)


 ### Creating inputs for models
  ETA <- list()
  
  # Model 1: Additive model (EL, without markers)
  ETA1 <- list(ENV=list(~factor(Y$ENV)-1,model='BRR'),
           VAR=list(~factor(Y$VAR)-1,model='BRR')
          )
  ETA$ETA1 <- ETA1
  
  # Model 2: Additive model (EG, with markers)
  myIndex <- match(Y$VAR,rownames(G))
  EVD <- eigen(G)
  PC <- EVD$vectors[,EVD$values>1e-5]
  for(i in 1:ncol(PC)) PC[,i] <- PC[,i]*sqrt(EVD$values[i])
  PC <- PC[myIndex,]
  ETA1$VAR <- list(X=PC,model='BRR')
  ETA$ETA2 <- ETA1

  # Model 3: Additive model (EGW, with markers and env. cov.)
  W=scale(W)/sqrt(ncol(W))
  ETA1$COV=list(X=W, model='BRR')
  ETA$ETA3 <- ETA1

  # Model 4: GxW Model 1 (EGW+GxW, Model 3 + interactions between markers and env. covariates)
  Omega <- tcrossprod(W)
  Omega <- Omega/mean(diag(Omega))
  GxW <- G[myIndex,myIndex]*Omega
  EVD <- eigen(GxW)
  PC <- EVD$vectors[,EVD$values>1e-5]
  for(i in 1:ncol(PC)) PC[,i] <- PC[,i]*sqrt(EVD$values[i])
  ETA1$GxW=list(X=PC, model='BRR')
  ETA$ETA4 <- ETA1
  
 ### Saving inputs
 save(X,Y,W,G,file="standardized_data.RData")
 save(ETA,file="ETA.RData")
  
```
[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
