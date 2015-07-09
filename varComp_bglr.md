### Variance Components & basic GxE analysis

  For the analyses involving markers and env. covariates we will use BGLR. Before que fit more complex models we compare the estimates of variance components of simpler models obtained with lme4 and BGLR.

   
```R
## Parameters
  genoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/X_2012_2014.rda'
  envCovFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/W_No_ctr_std.rda' 
  phenoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/Y.rda' 
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/varComp/'
 ###

 load(phenoFile)
 library(lme4)
 library(BGLR)

 ## lm4
 fm1=lmer(rdt~factor(YEAR)+(1|VAR)+(1|REGION/LOC)+(VAR|REGION),data=Y)  
 
 ## BGLR


 

```
|[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md):|
