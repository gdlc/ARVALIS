### Variance Components BGLR vs lmer

  For the analyses involving markers and env. covariates we will use BGLR. Before que fit more complex models we compare the estimates of variance components of simpler models obtained with lme4 and BGLR.

   
```R
## Parameters
  inputFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/standardized_data.RData'
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/varComp_bglr/'
 ###

 load(inputFile)
 library(lme4)
 library(BGLR)

 ## lm4
  Y$YEARxLOC=paste(Y$YEAR,"x",Y$LOC)
  Y$YEARxREGION=paste(Y$YEAR,"x",Y$REGION)
  
  
  fmLMER=lmer(rdt~(1|VAR)+(1|YEAR)+(1|REGION)+(1|LOC)+
               (1|YEARxREGION)+(1|YEARxLOC),data=Y)  
 
 
 ## BGLR
  LP=list(VAR=list(~factor(Y$VAR)-1,model='BRR') ,
          YEAR=list(~factor(Y$YEAR)-1,model='BRR'),
          REGION=list(~factor(Y$REGION)-1,model='BRR'),
          LOC=list(~factor(Y$LOC)-1,model='BRR'),
          YEARxREG=list(~factor(Y$YEARxREGION)-1,model='BRR'),
          YEARxLOC=list(~factor(Y$YEARxLOC)-1,model='BRR')
      )
      
      
  fmBGLR=BGLR(y=Y$rdt,ETA=LP,nIter=32000,burnIn=2000,thin=10)  
 
  summary(fmLMER)
  
  ## VARIANCE COMPONENTS BGLR
   fmBGLR$varE
   fmBGLR$ETA$VAR$varB
   fmBGLR$ETA$YEAR$varB
   fmBGLR$ETA$REGION$varB
   fmBGLR$ETA$LOC$varB
   fmBGLR$ETA$YEARxREG$varB
   fmBGLR$ETA$YEARxLOC$varB
   
  ## Trace plots
  
    plot(scan('ETA_VAR_varB.dat'),type='o',cex=.5,col=4)
    plot(scan('ETA_YEAR_varB.dat'),type='o',cex=.5,col=4) #....
    
    # Differences between lmer and BGLR are mostly due to the fact that BGLR reports posterior means not post. modes
    plot(density(scan('ETA_YEAR_varB.dat')),col=4) ; abline(v=fmBGLR$ETA$YEAR$varB,lty=2,col=2 )

    # Predictions are almost equivalent
    
    plot(predict(fmLMER),fmBGLR$yHat,cex=.5,col=4)
```
[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
