## Full-data  analysis
## Parameters
  genoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/X_2012_2014.rda'
  envCovFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/W_No_ctr_std.rda' 
  phenoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/Y.rda' 
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/full_data_models/'
 ###

 dir.create(outputFolder) 
 setwd(outputFolder)
 load(genoFile)
 load(envCovFile)
 OUT=matrix(nrow=5,ncol=5,NA)
 colnames(OUT)=c('Eror','G','E','GxE','LxE')


### Model 1: Additive model (without markers)

### Model 2: Additive model (with markers)

### Model 3: Additive model (with markers and env. cov.)

### Model 4: GxE Model 1 (Model 3 + interactions between markers and env. covariates).

### Model 5: GxE Model 2 (Model 3 + interactions between markers and environments).


[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
