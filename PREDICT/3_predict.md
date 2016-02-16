

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
 trnPC.G<-'/Users/gustavodeloscampos/ARVALIS/predict/output/PC_G.RData'
 trnPC.GW<-'/Users/gustavodeloscampos/ARVALIS/predict/output/PC_GW.RData'
 
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
 load(trnPC.G)
 load(trnPC.GW)
 load(tstEnvCov)
 load(tstGeno)
 load(tstPheno)

 
 ## Intercept
    yHat=rep(fm$mu+mean(fm$ETA$year$b),nrow(Y2)); var(yHat)
   
 ## Add location effect
   for(i in 1:length(yHat)){   
            tmp<-grep(x=names(fm$ETA$loc$b),pattern=Y2$LOC[i])
            tmp<-ifelse(length(tmp)==0,0,fm$ETA$loc$b[tmp])
   			yHat[i]=yHat[i]+tmp 
   	}
    var(yHat)
    
 ## Add genotype effect
   # Compute G21
   tmp<-which(colnames(X2)%in%colnames(X))
   X2=X2[,tmp]
   X2=scale(X2,scale=sdX,center=meansX)
   for(i in 1:ncol(X2)){
   	 tmp=is.na(X2[,i])
   	 if(any(tmp)){ X2[tmp,i]=0 }
   }
   G21<-tcrossprod(X2,X)/ncol(X)
   uHat=G21%*%GInv%*%PC.G%*%fm$ETA$G$b     
   for(i in 1:length(yHat)){   
       tmp<-grep(x=rownames(X2),pattern=Y2$VAR[i])
       tmp<-ifelse(length(tmp)==0,0,uHat[tmp])
   	   yHat[i]=yHat[i]+tmp 
   	}
   var(yHat)
   
 ## Add env. covariate effects
   tmp<-which(colnames(W2)%in%colnames(W))
   W2=W2[,tmp]
   W2=scale(W2,center=attributes(W)[[1]],scale=attributes(W)[[2]])
   W2<-W2/sqrt(ncol(W))
   tmp=W2%*%fm$ETA$W$b
   yHat=yHat+tmp   
   var(yHat)
   
   
 ## Add GW effect
  WW21<-tcrossprod(W2,W)
  ID.TST=as.integer(factor(x=Y2$VAR,levels=rownames(G21),ordered=T))
  ID.TRN= as.integer(factor(x=Y$VAR,levels=colnames(G21),ordered=T))
  GW21<-WW21*G21[ID.TST,ID.TRN]
  TMP=GW21%*%GWInv
  gwHat<-PC.GW%*%fm$ETA$GW$b
  tmp<-TMP%*%gwHat
  yHat=yHat+tmp
  
 ```
