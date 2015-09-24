### Variance Components 

In this module we use the lme4 package to fit a mixed effects model for yield. We fit three models models

   fm1:  y=YEAR(fixed) + VAR(random) + LOC(random) + Error 
   
   fm2:  y=YEAR(fixed) + VAR(random) + REGION(random) + LOC-within Region(random) + Error 
   
   fm3:  y=YEAR(fixed) + VAR(random) + REGION(random) + LOC-within Region(random) + YEARxREGION(random) + Error 
   
```R
 ## Parameters. 'inputFolder' was previously saved with X,Y,W and G matrices
  inputFolder='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/'
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/varComp_lmer/'
 ###
 dir.create(outputFolder)
 load(paste0(inputFolder,"/standardized_data.RData"))
 library(lme4)
 Y$YEARxLOC=paste(Y$YEAR,"x",Y$LOC)
 Y$YEARxREGION=paste(Y$YEAR,"x",Y$REGION)

  fm1 <- lmer(y~factor(YEAR)+(1|VAR)+(1|LOC),data=Y)  

  fm2 <- lmer(y~factor(YEAR)+(1|VAR)+(1|REGION/LOC),data=Y)  ## Partitioning varinace of LOC into Region and LOC within region.

  fm3 <- lmer(y~(1|VAR)+(1|YEAR)+(1|REGION)+(1|LOC)+
               (1|YEARxREGION)+(1|YEARxLOC),data=Y)  
  summary(fm1)
  summary(fm2)
  summary(fm3)
  save(fm3,file=paste0(outputFolder,'fm.RData'))
  
  colMeans(X)
  
```
[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
