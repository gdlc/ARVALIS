## Full-data  analysis

```R
 ## Parameters
  genoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/X_2012_2014.rda'
  envCovFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/W_No_ctr_std.rda' 
  phenoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/Y.rda' 
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/full_data_models/'
 
  library(BGLR)
  nIter=2000; burnIn=500
 ###
 
 dir.create(outputFolder) 
 setwd(outputFolder)
 load(genoFile) ; load(envCovFile); load(phenoFile)
 OUT=matrix(nrow=5,ncol=4,NA)
 colnames(OUT)=c('Eror','G','E','GxE')

 ##
   IDs=intersect(Y$VAR,rownames(X))
   Y=Y[Y$VAR%in%IDs,]
   W=W_No_ctr_std[Y$VAR%in%IDs,]
   X=X[IDs,];   stopifnot(all(Y$VAR%in%rownames(X))); stopifnot(all(rownames(X)%in%Y$VAR)) 

 ### Model 1: Additive model (without markers)
  ETA=list(ENV=list(~factor(Y$ENV)-1,model='BRR'),
           VAR=list(~factor(Y$VAR)-1,model='BRR')
          )
  fmL=BGLR(y=Y$rdt,ETA=ETA,saveAt='EL_',nIter=nIter,burnIn=burnIn)
 
 
 ### Model 2: Additive model (with markers)
  G <- tcrossprod(scale(X))/ncol(X)
  myIndex=as.integer(factor(x=Y$VAR,ordered=T,levels=rownames(G)))
  myIndex <- match(Y$VAR,rownames(G))
  SVD <- eigen(G)
  PC <- SVD$vectors
  for(i in 1:ncol(PC)) PC[,i] <- PC[,i]*sqrt(SVD$values)
  PC=PC[myIndex,]
  ETA$VAR=list(X=PC,model='BRR')
  fmLM=BGLR(y=Y$rdt,ETA=ETA,saveAt='EM_',nIter=nIter,burnIn=burnIn)

 
 ### Model 3: Additive model (with markers and env. cov.)
  W=scale(W)
  ETA$COV=list(X=W, model='BRR')
  
  
 ### Model 4: GxE Model 1 (Model 3 + interactions between markers and env. covariates).

 ### Model 5: GxE Model 2 (Model 3 + interactions between markers and environments).

```
 


[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
