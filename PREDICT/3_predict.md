

## Module 2: Model fitting


```R

rm(list=ls())
 library(BGLR)
 library(BGData)
 
## Change parameters here
 outputFolder='/Users/gustavodeloscampos/ARVALIS/predict/output'
 fittedModel='/Users/gustavodeloscampos/ARVALIS/predict/output/fm.RData'
 trnEnvCov<-'/Users/gustavodeloscampos/ARVALIS/predict/output/W.RData'
 trnGeno<-'/Users/gustavodeloscampos/ARVALIS/predict/output/X.RData'
 trnPheno<-'/Users/gustavodeloscampos/ARVALIS/predict/output/Y.RData'
 trnEVDG<-'/Users/gustavodeloscampos/ARVALIS/predict/output/EVD_G.RData'
 trnPCGW<-
 
 tstEnvCov<-'/Users/gustavodeloscampos/ARVALIS/predict/newData/W2.RData'
 tstGeno<-'/Users/gustavodeloscampos/ARVALIS/predict/newData/X2.RData'
 tstPheno<-'/Users/gustavodeloscampos/ARVALIS/predict/newData/Y2.RData'
 
 
 trnScalesAndMeans<-'/Users/gustavodeloscampos/ARVALIS/predict/output/meansAndSDs.RData'
 trnGInv<-'/Users/gustavodeloscampos/ARVALIS/predict/output/G_GInv.RData'
 trnGWInv<-'/Users/gustavodeloscampos/ARVALIS/predict/output/GW_GWInv.RData'
###
 
 
 dir.create(outputFolder)
 setwd(outputFolder)
 
 setwd(outputFolder)
 load(fittedModel)
 load(trnEnvCov)
 load(trnGeno)
 load(trnPheno)
 load(trnScalesAndMeans)
 load(trnGInv); rm(G)
 load(trnGWInv); rm(GW)
 load(tstEnvCov)
 load(tstGeno)
 load(tstPheno)

 
 ## Intercept
   yHat=rep(fm$mu,nrow(Y2))
 ## Add location effect
   for(i in 1:length(yHat)){   yHat=yHat+fm$loc$b[grep(x=names(fm$loc$b),pattern=Y$LOC)] }
   
 ## Add genotype effect
   # Compute G21
   tmp<-which(colnames(X2)%in%colnames(X))
   X2=X2[,tmp]
   G21<-tcrossprod(X2,X)/sqrt(ncol(X))
   uHat=G21%*%GInv%*%PC.G%*%fm$G$b
 
 ## Add env. covariate effects
 
 ## Add GW effect
 
 ```
