### Structure of Genotypes and of Env. Covariates Based on Eigen-decomposition


```R
## Parameters
  genoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/X_2012_2014.rda'
  envCovFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/W_No_ctr_std.rda' 
  phenoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/Y.rda' 
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/eigen/'
 ###
 
 dir.create(outputFolder) 
 setwd(outputFolder)
 load(genoFile)
 load(envCovFile)
 
 pdf('eigenAnalysis.pdf')
 
 ### Gentoypes ###################################
 
  # MAF and missing values frequency
   freqNA <- apply(X,2,function(x) sum(!x%in%c(0,1,2))/length(x)) 
   MAF <- apply(X,2,function(x) mean(x[x%in%c(0,1,2)])/2)
   MAF <- ifelse(MAF>0.5,1-MAF,MAF)

  # Remove markers with MAF < 0.05 and freqNA > 0.5
   index <- which(MAF<0.05 | freqNA>0.5)
   if(length(index)>0) X <- X[,-index]
 
  # Genomic relationship matrix
   Z <- scale(X,center=TRUE,scale=TRUE)
   G <- tcrossprod(Z)/ncol(Z)
 
  # Principal component analysis
   EVD<- eigen(G)
   varExp <- 100*EVD$values/sum(EVD$values)
   par(mfrow=c(1,2))
   plot(EVD$vectors[,1:2],main="Genotypes: PC 1 vs PC 2", xlab="PC 1",ylab="PC 2",col=4,cex=.5)
   plot(cumsum(varExp[1:20]),main="Genotypes: cumulative variance explained",
        ylab="% of Total",xlab="PC (Top 20)",type='o',col=4,cex=.5)
   abline(h=60,lty=2,col="red")
  
  ## Environmental covariates ################
  
  
 dev.off()
```
[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)

[See PDF PrinComp] (https://github.com/gdlc/ARVALIS/blob/master/PC_analysis_G.pdf)
