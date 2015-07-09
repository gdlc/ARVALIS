### Eigen analyses of genotypes and of env. covariates.

   
```R
## Parameters
  genoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/X_2012_2014.rda'
  envCovFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/W_No_ctr_std.rda' 
  phenoFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/getData/output/Y.rda' 
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/eigen/'
 ###
 
 dir.create(outputFolder) 
 load(genoFile)
 load(envCovFile)
 
 # MAF and missing values frequency
 freqNA <- apply(X,2,function(x) sum(!x%in%c(0,1,2))/length(x)) 
 MAF <- apply(X,2,function(x) mean(x[x%in%c(0,1,2)])/2)
 MAF <- ifelse(MAF>0.5,1-MAF,MAF)

 # Remove markers with MAF < 0.05 and freqNA > 0.5
 index <- which(MAF<0.05 | freqNA>0.5)
 if(length(index)>0) X <- X[,-index]
 
 # Genomic matrix
 Z <- scale(X,center=TRUE,scale=TRUE)
 G <- tcrossprod(Z)/ncol(Z)
 
 # Principal component analysis
 PC <- princomp(G)
 varExp <- PC$sdev^2
 varExp <- 100*varExp/sum(varExp)
 par(mfrow=c(1,2))
 plot(PC$scores[,1],PC$scores[,2],main="PC 1 vs PC 2", xlab="PC 1",ylab="PC 2")
 plot(cumsum(varExp[1:20]),main="Cumulated explained variance",ylab="% of Total",xlab="PC (Top 20)")
 abline(h=60,lty=2,col="red")

```
[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
[See PDF PrinComp] (https://github.com/gdlc/ARVALIS/blob/master/PC_analysis_G.pdf)
