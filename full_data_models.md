## Full-data  analysis

```R
 ## Parameters
  genoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/X_2012_2014.rda'
  envCovFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/W_No_ctr_std.rda' 
  phenoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/Y.rda' 
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/full_data_models/'
 
  library(BGLR)
  nIter=1200; burnIn=200
 ###
 
 dir.create(outputFolder) 
 setwd(outputFolder)
 load(genoFile) ; load(envCovFile); load(phenoFile)
 OUT=matrix(nrow=4,ncol=5,NA)
 colnames(OUT)=c('E','G','W','GxW','Error')
 rownames(OUT)=c('EL','EG','EGW','EGW_GxW')
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
  fmEL=BGLR(y=Y$rdt,ETA=ETA,saveAt='EL_',nIter=nIter,burnIn=burnIn)
  OUT['EL','E']=fmEL$ETA$ENV$varB
  OUT['EL','G']=fmEL$ETA$VAR$varB
  OUT['EL','Error']=fmEL$varE
 
 ### Model 2: Additive model (with markers)
  G <- tcrossprod(scale(X))
  G<-G/mean(diag(G))
  myIndex=as.integer(factor(x=Y$VAR,ordered=T,levels=rownames(G)))
  myIndex <- match(Y$VAR,rownames(G))
  SVD <- eigen(G)
  PC <- SVD$vectors
  for(i in 1:ncol(PC)) PC[,i] <- PC[,i]*sqrt(SVD$values)
  PC=PC[myIndex,]
  ETA$VAR=list(X=PC,model='BRR')
  fmEG=BGLR(y=Y$rdt,ETA=ETA,saveAt='EG_',nIter=nIter,burnIn=burnIn)
  OUT['EG','E']=fmEG$ETA$ENV$varB
  OUT['EG','G']=fmEG$ETA$VAR$varB
  OUT['EG','Error']=fmEG$varE

 ### Model 3: Additive model (with markers and env. cov.)
  W=scale(W)/sqrt(ncol(W))
  ETA$COV=list(X=W, model='BRR')
  fmEGW=BGLR(y=Y$rdt,ETA=ETA,saveAt='EGW_',nIter=nIter,burnIn=burnIn)
  OUT['EGW','E']=fmEGW$ETA$ENV$varB
  OUT['EGW','G']=fmEGW$ETA$VAR$varB
  OUT['EGW','W']=fmEGW$ETA$COV$varB
  OUT['EGW','Error']=fmEGW$varE
  
 ### Model 4: GxW Model 1 (Model 3 + interactions between markers and env. covariates)
  Omega <- tcrossprod(W)
  Omega<-Omega/mean(diag(Omega))
  GxW <- G[myIndex,myIndex]*Omega
  EVD <- eigen(GxW)
  PC<-EVD$vectors[,EVD$values>1e-5]
  for(i in 1:ncol(PC)){
    PC[,i]=EVD$vectors[,i]*sqrt(EVD$values[i])
  }
  
  ETA$GxW=list(X=PC, model='BRR')
  fmEGW_GxW=BGLR(y=Y$rdt,ETA=ETA,saveAt='EGW_GxW_',nIter=nIter,burnIn=burnIn)
  OUT['EGW_GxW','E']=fmEGW_GxW$ETA$ENV$varB
  OUT['EGW_GxW','G']=fmEGW_GxW$ETA$VAR$varB
  OUT['EGW_GxW','W']=fmEGW_GxW$ETA$COV$varB
  OUT['EGW_GxW','GxW']=fmEGW_GxW$ETA$GxW$varB
  OUT['EGW_GxW','Error']=fmEGW_GxW$varE
  
  round(OUT,3)

```
 
Estimated Variace Components


| Model     | E    |  G     | W    | GxW   | Error |
| --------- |-----:| -----:| -----:| -----:| -----:| 
EL      | 201.562  | 14.624  |     --  |    --  | 22.651 | 
EG      |  201.738 |  24.751  |     --  |    --  | 22.701 | 
EGW      | 143.213 |  25.093 |  40.783 |     --  | 22.628 | 
EGW_GxW  | 146.694  | 21.254 |  34.464 |  6.125 |  18.305 | 


[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
