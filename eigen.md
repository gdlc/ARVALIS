### Structure of Genotypes and of Env. Covariates Based on Eigen-decomposition

```SAS
 proc reg;
 model y=x1+x2;
 run;
```

```R
## Parameters
  inputFolder='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/'
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/eigen/'
 ###
 
 dir.create(outputFolder) 
 setwd(outputFolder)
 load(paste0(inputFolder,"/standardized_data.RData"))
 
 pdf('eigenAnalysis.pdf')
 ### Genotypes ###################################

  # Principal component analysis
   EVD <- eigen(G)
   
  # Eigenvalues
   varExp <- 100*EVD$values/sum(EVD$values)
   plot(x=0:20,y=cumsum(c(0,varExp[1:20])),main="Genotypes: cumulative variance explained",
        ylab="% of Total",xlab="PC (Top 20)",type='o',col=4,cex=.5)
   abline(h=50,lty=2,col="red")
   abline(a=0,b=1,col=8,lty=2)

   # Eigenvectors
   plot(EVD$vectors[,1:2],main="Genotypes: PC 1 vs PC 2", xlab="PC 1",ylab="PC 2",col=4,cex=.5)
   
  ## Environmental covariates ################
   SVD <- svd(W,nu=nrow(W),nv=ncol(W)-1)
   varExp <- 100*SVD$d/sum(SVD$d)

  # Eigenvalues
   plot(x=0:20,c(0,cumsum(varExp[1:20])),main="Env. Covariates: cumulative variance explained",
        ylab="% of Total",xlab="PC (Top 20)",type='o',col=4,cex=.5)
   abline(h=60,lty=2,col="red")

  # Eigenvectors
   plot(SVD$u[,1:2],main="Env. Covariates: PC 1 vs PC 2 (by location)", xlab="PC 1",ylab="PC 2",col="white",cex=.5)
   labels <- unique(Y$LOC)
   for(i in 1:length(labels)){
      tmp=Y$LOC==labels[i]
      points(x=SVD$u[tmp,1],y=SVD$u[tmp,2],col=i,cex=.5)
   }

  plot(SVD$u[,1:2],main="Env. Covariates: PC 1 vs PC 2 (by region)", xlab="PC 1",ylab="PC 2",col="white",cex=.5)
   labels <- unique(Y$REGION)
   for(i in 1:length(labels)){
      tmp=Y$REGION==labels[i]
      points(x=SVD$u[tmp,1],y=SVD$u[tmp,2],col=i,cex=.5)
   }
   
 dev.off()
```
[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)

Plots produced by the above-script.

![ScreenShot](https://github.com/gdlc/ARVALIS/blob/master/eigenAnalysis.pdf)
