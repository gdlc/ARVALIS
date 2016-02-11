
##
# This module prepare data for: 
#      (i) fitting the model to a training data set
#     (ii) deriving predictions for testing data
##


rm(list=ls())
 library(BGLR)
 library(BGData)
 
## Change parameters here
 phenotypeFile<-'~/Dropbox/ARVALIS/data/Data2015/Y.RData'
 envCovFile<-'~/Dropbox/ARVALIS/data/Data2015/W.RData'
 genFile<-'~/Dropbox/ARVALIS/data/Data2015/X.RData'
 outputFolder='/Users/gustavodeloscampos/ARVALIS/predict/output'
 nMinLines<-20
 nMinEnv<-20
###
 
 
 dir.create(outputFolder)
 setwd(outputFolder)
 load(phenotypeFile)
 load(envCovFile)
 load(genFile)
 
 
 for(i in 1:5){
	 tmp<-names(table(Y$VAR))[table(Y$VAR)>=nMinLines]
 	W<-W[Y$VAR%in%tmp,]
 	Y<-Y[Y$VAR%in%tmp,]
 
 
 	tmp<-names(table(Y$ENV))[table(Y$ENV)>=nMinEnv]
 	W<-W[Y$ENV%in%tmp,]
 	Y<-Y[Y$ENV%in%tmp,]
 
 
 	tmp=rownames(X)%in%Y$VAR
 	X<-X[tmp,]
 }
 
  
 ## Removing markers with more than 5% of NAs or MAF<3%
 
  tmp=apply(FUN=function(x) mean(is.na(x)), MARGIN=2,X=X)
  X<-X[,tmp<.05]
  tmp=colMeans(X,na.rm=T)/2
  maf<-ifelse(tmp>.5,1-tmp,tmp)
  X<-X[,maf>.03]
 
 
 ## Centers and scales are saved to be used in predictions
  sdX<-apply(X=X,MARGIN=2,FUN=sd,na.rm=T)*sqrt((nrow(X)-1)/nrow(X))
  meansX<-colMeans(X,na.rm=T)
  meansW<-colMeans(W,na.rm=T)
  sdW<-apply(X=W,MARGIN=2,FUN=sd,na.rm=T)*sqrt((nrow(W)-1)/nrow(W))
  save(sdX,sdW,meansW,meansX,file='meansAndSDs.RData')


  # centering and scaling markers and covariates
  W=scale(W,scale=sdW,center=meansW)
  W<-W/sqrt(ncol(W))
  

  
  for(i in 1:ncol(X)){
	xi=X[,i]
	tmp=is.na(xi)
	xi=(xi-meansX[i])/sdX[i]
	if(any(tmp)){ xi[tmp]=0 }
	X[,i]=xi
	print(i)
  }
  
  save(X,file='X.RData')
  save(W,file='W.RData')
  save(Y,file='Y.RData')
  
  G<-tcrossprod(X)
  G<-G/mean(diag(G))

  EVD.G<-eigen(G)
  EVD.G$vectors=EVD.G$vectors[,EVD.G$values>1e-5]
  EVD.G$values=EVD.G$values[EVD.G$values>1e-5]
  save(EVD.G,file='EVD_G.RData')
  
  
  ## continue here
  PC.G<-EVD.G$vectors
  rownames(PC.G)<-rownames(G)
  for(i in 1:ncol(PC.G)){ PC.G[,i]<-PC.G[,i]*sqrt(EVD.G$values[i]) }
  
  
  TMP<-EVD.G$vectors; for(i in 1:ncol(EVD.G$vectors)){ TMP[,i]=TMP[,i]/sqrt(EVD.G$values[i]) }
  
  GInv<-tcrossprod(TMP)
   
  save(G,GInv,file='G_GInv.RData')

  WW<-tcrossprod(W)
 
  ID.VAR.NUM<-as.integer(factor(x=Y$VAR,ordered=T,levels=rownames(G)))
  stopifnot( all.equal(rownames(G)[ID.VAR.NUM],Y$VAR))
  
  ZGZ<-G[ID.VAR.NUM,ID.VAR.NUM]
  GW<-WW*ZGZ

  EVD.GW<-eigen(GW)
  EVD.GW$vectors=EVD.GW$vectors[,EVD.GW$values>1e-5]
  EVD.GW$values=EVD.GW$values[EVD.GW$values>1e-5]

  PC.GW<-EVD.GW$vectors
  for(i in 1:ncol(PC.GW)){ PC.GW[,i]<-PC.GW[,i]*sqrt(EVD.GW$values[i]) }
  save(PC.GW,file='PC_GW.RData')
  
  TMP<-EVD.GW$vectors; for(i in 1:ncol(EVD.GW$vectors)){ TMP[,i]=TMP[,i]/sqrt(EVD.GW$values[i]) }
  
  GWInv<-tcrossprod(TMP)
  
  save(GW,GWInv,file='GW_GWInv.RData')

  
  
  
  
  
  
  
  
