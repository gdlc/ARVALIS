## Cross validation 1 and 2

**Cross validation 1**.
Newly developed lines. Random assigns lines to folds. All the records for a given line are assigned to the same fold.
 
**Cross validation 2**.
To assess the ability of models to predict the performance of lines using the data collected in other environments.
 
```R
 ## Parameters
  inputFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/standardized_data.RData'
  inputETA='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/ETA.RData'
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/cross_validation/'
  CV <- 1
  library(BGLR)
  nIter=1200; burnIn=200
 ###
 
 dir.create(outputFolder) 
 setwd(outputFolder)
 load(inputFile)
 
 nfolds <- 5
 IDs <- rownames(X)
 IDy <- as.character(Y$VAR)
 folds <- rep(1:nfolds,each=ceiling(length(IDs)/nfolds))[order(runif(length(IDs)))]
 CVfolds <- rep(NA,nrow(Y))

 if(CV==1)
 {
   for(i in 1:nfolds)
   {
       tmp <- IDy%in%IDs[fold==i]
       CVfolds[tmp] <- i
   }
 }

 if(CV==2)
 {
   for(i in IDs)
   {
       tmp <- which(IDy==i)
       ni <- length(tmp)
       fold0 <- sample(1:nfolds,size=ni,replace=length(tmp)>nfolds)
       CVfolds[tmp] <- fold0
   }
 }
 
 ### Model 1: Additive model (without markers)
  load(inputETA)
  ETA <- ETA$ETA1
  fmEL <- BGLR(y=Y$rdt,ETA=ETA,saveAt='EL_',nIter=nIter,burnIn=burnIn)
  OUT['EL','E']=fmEL$ETA$ENV$varB
  OUT['EL','G']=fmEL$ETA$VAR$varB
  OUT['EL','Error']=fmEL$varE
 
 ### Model 2: Additive model (with markers)
  load(inputETA)
  ETA <- ETA$ETA2
  fmEG <- BGLR(y=Y$rdt,ETA=ETA,saveAt='EG_',nIter=nIter,burnIn=burnIn)
  OUT['EG','E']=fmEG$ETA$ENV$varB
  OUT['EG','G']=fmEG$ETA$VAR$varB
  OUT['EG','Error']=fmEG$varE

 ### Model 3: Additive model (with markers and env. cov.)
  load(inputETA)
  ETA <- ETA$ETA3
  fmEGW <- BGLR(y=Y$rdt,ETA=ETA,saveAt='EGW_',nIter=nIter,burnIn=burnIn)
  OUT['EGW','E']=fmEGW$ETA$ENV$varB
  OUT['EGW','G']=fmEGW$ETA$VAR$varB
  OUT['EGW','W']=fmEGW$ETA$COV$varB
  OUT['EGW','Error']=fmEGW$varE
  
 ### Model 4: GxW Model 1 (Model 3 + interactions between markers and env. covariates)
  load(inputETA)
  ETA <- ETA$ETA4
  fmEGW_GxW <- BGLR(y=Y$rdt,ETA=ETA,saveAt='EGW_GxW_',nIter=nIter,burnIn=burnIn)
  OUT['EGW_GxW','E']=fmEGW_GxW$ETA$ENV$varB
  OUT['EGW_GxW','G']=fmEGW_GxW$ETA$VAR$varB
  OUT['EGW_GxW','W']=fmEGW_GxW$ETA$COV$varB
  OUT['EGW_GxW','GxW']=fmEGW_GxW$ETA$GxW$varB
  OUT['EGW_GxW','Error']=fmEGW_GxW$varE
  
  round(OUT,3)
```
