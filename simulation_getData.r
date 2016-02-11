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
 
 tmp<-names(table(Y$VAR))[table(Y$VAR)>=nMinLines]
 W<-W[Y$VAR%in%tmp,]
 Y<-Y[Y$VAR%in%tmp,]
 
 
 tmp<-names(table(Y$ENV))[table(Y$ENV)>=nMinEnv]
 W<-W[Y$ENV%in%tmp,]
 Y<-Y[Y$ENV%in%tmp,]
 
 
 tmp=rownames(X)%in%Y$VAR
 X<-X[tmp,]
 
 ## Removing markers with too many NAs or MAF<3%
 
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
  G<-tcrossprod(X)
  G<-G/mean(diag(G))

  EVD.G<-eigen(G)
  EVD.G$vectors=EVD.G$vectors[,EVD.G$values>1e-5]
  EVD.G$values=EVD.G$values[EVD.G$values>1e-5]

  ## continue here
  PC.G<-EVD.G$vectors
  for(i in 1:ncol(PC.G)){ PC.G[,i]<-PC.G[,i]*sqrt(EVD.G$values[i]) }

  TMP<-EVD.G$vectors; for(i in 1:ncol(EVD.G$vectors)){ TMP[,i]=TMP[,i]/sqrt(EVD.G$values[i]) }
  
  GInv<-tcrossprod(TMP)
   
  save(G,GInv,file='G_GInv.RData')

  O<-tcrossprod(W)
 
  ID.VAR.NUM<-as.integer(factor(x=Y$VAR,ordered=T,levels=rownames(G)))
  stopifnot( all.equal(rownames(G)[ID.VAR.NUM],Y$VAR))
  
  ZGZ<-G[ID.VAR.NUM,ID.VAR.NUM]
  OG<-O*ZGZ

  EVD.OG<-eigen(OG)
  EVD.OG$vectors=EVD.OG$vectors[,EVD.OG$values>1e-5]
  EVD.OG$values=EVD.OG$values[EVD.OG$values>1e-5]

  ## continue here
  PC.OG<-EVD.OG$vectors
  for(i in 1:ncol(PC.OG)){ PC.OG[,i]<-PC.OG[,i]*sqrt(EVD.OG$values[i]) }

  TMP<-EVD.OG$vectors; for(i in 1:ncol(EVD.OG$vectors)){ TMP[,i]=TMP[,i]/sqrt(EVD.OG$values[i]) }
  
  OGInv<-tcrossprod(TMP)
  
  save(OG,OGInv,file='OG_OGInv.RData')

  
  
  
  
  
  
  
  