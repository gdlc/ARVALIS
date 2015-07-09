### Variance Components 

In this module we use the lme4 package to fit a mixed effects model for yield. We fit three models models

   fm1:  y=YEAR(fixed) + VAR(random) + LOC(random) + Error 
   
   fm2:  y=YEAR(fixed) + VAR(random) + REGION(random) + LOC-within Region(random) + Error 
   
   fm3:  y=YEAR(fixed) + VAR(random) + REGION(random) + LOC-within Region(random) + YEARxREGION(random) + Error 
   
```R
## Parameters
  genoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/X_2012_2014.rda'
  envCovFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/W_No_ctr_std.rda' 
  phenoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/Y.rda' 
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/varComp_lmer/'
 ###
 dir.create(outputFolder)
 load(phenoFile)
 library(lme4)
 Y$YEARxLOC=paste(Y$YEAR,"x",Y$LOC)
 Y$YEARxREGION=paste(Y$YEAR,"x",Y$REGION)

  fm1=lmer(rdt~factor(YEAR)+(1|VAR)+(1|LOC),data=Y)  

  fm2=lmer(rdt~factor(YEAR)+(1|VAR)+(1|REGION/LOC),data=Y)  ## Partitioning varinace of LOC into Region and LOC within region.

  fm3=lmer(rdt~(1|VAR)+(1|YEAR)+(1|REGION)+(1|LOC)+
               (1|YEARxREGION)+(1|YEARxLOC),data=Y)  
  summary(fm1)
  summary(fm2)
  summary(fm3)
  save(fm3,file=paste0(outputFolder,'fm.RData'))
  
```
[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
