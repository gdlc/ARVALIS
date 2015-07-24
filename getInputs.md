This code is for preparing data: matching genotypic and genotypic information, applying a quality control to genotypes.

```R
## Parameters
  genoFile='/Users/epidemiology/Documents/Biostatistics/ARVALIS/PIPELINE_2014/data/X_2012_2014.rda'
  envCovFile='/Users/epidemiology/Documents/Biostatistics/ARVALIS/PIPELINE_2014/data/W_No_ctr_std.rda' 
  phenoFile='/Users/epidemiology/Documents/Biostatistics/ARVALIS/PIPELINE_2014/data/Y.rda' 
  outputFolder='/Users/epidemiology/Documents/Biostatistics/ARVALIS/PIPELINE_2014/input/'
  
  thresholdMAF=0.05
  thresholdNAfreq=0.3
  valid_markers = c(0,1,2)
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
```
[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
