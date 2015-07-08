###Estimates variance components using lme4

```R
## Parameters
  genoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/X_2012_2014.rda'
  envCovFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/W_No_ctr_std.rda' 
  phenoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/Y.rda' 
 ###


 load(phenoFile)
  
 library(lme4)
  
 ## Mixed model with Year as fixed effect, Var, REGION and LOC-within-REGION as random
  
  fm=lmer(rdt~factor(YEAR)+(1|VAR)+(1|REGION/LOC),data=Y)
  summary(fm)
```R
