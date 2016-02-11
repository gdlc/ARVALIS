

## Module 2: Model fitting


```R


rm(list=ls())
 library(BGLR)
 library(BGData)
 
## Change parameters here
 nIter=80000
 burnIn=20000
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
  Z.LOC=as.matrix(model.matrix(~factor(Y$YEAR))[,-1])
  Z.YEARxLOC=as.matrix(model.matrix(~factor(paste(Y$YEAR,Y$LOC,sep='-')))[,-1])
  
  ETA=list(
  			year=list(X=Z.YEAR,model='FIXED'),
  			loc=list(X=Z.LOC,model='BRR',saveEffects=T),
  			yearLoc=list(X=Z.YEARxLOC,model='BRR',saveEffects=T),
            G=list(X=ZPC.G,model='BRR',saveEffects=T),
            W=list(X=W,model='BRR',saveEffects=T),
            GW=list(X=PC.GW,model='BRR',saveEffects=T)
  		  )
  fm=BGLR(y=Y$rdt,ETA=ETA,nIter=nIter,burnIn=burnIn)
  save(fm,file='fm.RData')

```R

** Post-hoc ***

```R

  BYear=as.matrix(read.table('ETA_year_b.dat',header=T))
  BLoc=readBinMat('ETA_loc_b.bin')
  BYearLoc=readBinMat('ETA_yearLoc_b.bin')
  BG=readBinMat('ETA_G_b.bin')
  BW=readBinMat('ETA_W_b.bin')
  BGW=readBinMat('ETA_GW_b.bin')
   
  yHatYear=Z.YEAR%*%t(BYear)
   plot(apply(X=yHatYear,MARGIN=2,FUN=var))
  
  yHatLoc=Z.LOC%*%t(BLoc)
   plot(apply(X=yHatLoc,MARGIN=2,FUN=var))

  yHatYearLoc=Z.YEARxLOC%*%t(BYearLoc)
   plot(apply(X=yHatYearLoc,MARGIN=2,FUN=var))
  
  yHatG=ZPC.G%*%t(BG)
   plot(apply(X=yHatG,MARGIN=2,FUN=var))
  
  yHatW=W%*%t(BW)
   plot(apply(X=yHatW,MARGIN=2,FUN=var))
  
  yHatGW<-tcrossprod(PC.GW,BGW)
  
    plot(apply(X=yHatGW,MARGIN=2,FUN=var))
```
