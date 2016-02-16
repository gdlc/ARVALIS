

## Module 2: Model fitting


```R


rm(list=ls())
 library(BGLR)
 library(BGData)
 
## Change parameters here
 nIter=230000
 burnIn=30000
 outputFolder='/Users/gustavodeloscampos/ARVALIS/predict/output'
 phenotypeFile<-'/Users/gustavodeloscampos/ARVALIS/predict/output/Y.RData'
 envCovFile<-'/Users/gustavodeloscampos/ARVALIS/predict/output/W.RData'
 PC_GFile<-'/Users/gustavodeloscampos/ARVALIS/predict/output/PC_G.RData'
 PC_GWFile<-'/Users/gustavodeloscampos/ARVALIS/predict/output/PC_GW.RData'
 

###
 
 
 dir.create(outputFolder)
 setwd(outputFolder)
 load(phenotypeFile)
 load(envCovFile)
 load( PC_GFile)
 load( PC_GWFile)


  
  ID.VAR.NUM<-as.integer(factor(x=Y$VAR,ordered=T,levels=rownames(PC.G)))
  stopifnot( all.equal(rownames(G)[ID.VAR.NUM],Y$VAR))
  
  ZPC.G<-PC.G[ID.VAR.NUM,]
  
  Z.YEAR=as.matrix(model.matrix(~factor(Y$YEAR))[,-1])
  Z.LOC=as.matrix(model.matrix(~factor(Y$LOC))[,-1])
  Z.YEARxLOC=as.matrix(model.matrix(~factor(paste(Y$YEAR,Y$LOC,sep='-')))[,-1])
  
  ETA=list(
  			year=list(X=Z.YEAR,model='BRR'),
  			loc=list(X=Z.LOC,model='BRR'),
  			yearLoc=list(X=Z.YEARxLOC,model='BRR'),
            G=list(X=ZPC.G,model='BRR'),
            W=list(X=W,model='BRR'),
            GW=list(X=PC.GW,model='BRR')
  		  )
  fm=BGLR(y=Y$rdt,ETA=ETA,nIter=nIter,burnIn=burnIn)
  save(fm,file='fm.RData')



```R
