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
   index <- Y$VAR%in%IDs
   Y=Y[index,]
   W=W_No_ctr_std[index,]
   X=X[IDs,];   stopifnot(all(Y$VAR%in%rownames(X))); stopifnot(all(rownames(X)%in%Y$VAR)) 

 ### Model 1: Additive model (without markers)
  ETA=list(ENV=list(~factor(Y$ENV)-1,model='BRR'),
           VAR=list(~factor(Y$VAR)-1,model='BRR')
          )
  fm1=BGLR(y=Y$rdt,ETA=ETA,saveAt='EL_',nIter=nIter,burnIn=burnIn)
 
 
 ### Model 2: Additive model (with markers)
  G <- tcrossprod(scale(X))/ncol(X)
  myIndex=as.integer(factor(x=Y$VAR,ordered=T,levels=rownames(G)))
  myIndex <- match(Y$VAR,rownames(G))
  SVD <- eigen(G)
  PC <- SVD$vectors
  for(i in 1:ncol(PC)) PC[,i] <- PC[,i]*sqrt(SVD$values)
  PC=PC[myIndex,]
  ETA$VAR=list(X=PC,model='BRR')
  fm2=BGLR(y=Y$rdt,ETA=ETA,saveAt='EM_',nIter=nIter,burnIn=burnIn)

 
 ### Model 3: Additive model (with markers and env. cov.)
  W=scale(W)
  ETA$COV=list(X=W, model='BRR')
  fm3=BGLR(y=Y$rdt,ETA=ETA,saveAt='EMW_',nIter=nIter,burnIn=burnIn)
  
  
 ### Model 4: GxE Model 1 (Model 3 + interactions between markers and env. covariates)
  S <- tcrossprod(W)
  GW <- G[myIndex,myIndex]*S
  SVD <- eigen(GW)
  ETA$WG=list(V=SVD$vectors,d=SVD$values, model='RKHS')
  fm4=BGLR(y=Y$rdt,ETA=ETA,saveAt='EMW1_',nIter=nIter,burnIn=burnIn)


 ### Model 5: GxE Model 2 (Model 3 + interactions between markers and environments)


```
 


[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
