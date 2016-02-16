

## Module 2: Model fitting


```R

rm(list=ls())

timeIn<-proc.time()

 library(BGLR)
 library(BGData)
 
 outputFolder='/Users/gustavodeloscampos/ARVALIS/predict/output'
 fittedModel='/Users/gustavodeloscampos/ARVALIS/predict/output/fm.RData'
 trnEnvCov<-'/Users/gustavodeloscampos/ARVALIS/predict/output/W.RData'
 trnGeno<-'/Users/gustavodeloscampos/ARVALIS/predict/output/X.RData'
 trnPheno<-'/Users/gustavodeloscampos/ARVALIS/predict/output/Y.RData'
 trnEVDG<-'/Users/gustavodeloscampos/ARVALIS/predict/output/EVD_G.RData'
 trnPC.G<-'/Users/gustavodeloscampos/ARVALIS/predict/output/PC_G.RData'
 trnPC.GW<-'/Users/gustavodeloscampos/ARVALIS/predict/output/PC_GW.RData'
 
 tstEnvCov<-'/Users/gustavodeloscampos/ARVALIS/predict/newData/W2.RData'
 tstGeno<-'/Users/gustavodeloscampos/ARVALIS/predict/newData/X2.RData'
 tstPheno<-'/Users/gustavodeloscampos/ARVALIS/predict/newData/Y2.RData'
 
 
 trnScalesAndMeans<-'/Users/gustavodeloscampos/ARVALIS/predict/output/meansAndSDs.RData'
 trnGInv<-'/Users/gustavodeloscampos/ARVALIS/predict/output/G_GInv.RData'
 trnGWInv<-'/Users/gustavodeloscampos/ARVALIS/predict/output/GW_GWInv.RData'
###
 
 
 dir.create(outputFolder)
 setwd(outputFolder)
 
 setwd(outputFolder)
 load(fittedModel)
 load(trnEnvCov)
 load(trnGeno)
 load(trnPheno)
 load(trnScalesAndMeans)
 load(trnGInv); rm(G)
 load(trnGWInv); rm(GW)
 load(trnPC.G)
 load(trnPC.GW)
 load(tstEnvCov)
 load(tstGeno)
 load(tstPheno)


  nTST<-nrow(Y2)
 
 ## Intercept
    muHat=rep(fm$mu+mean(fm$ETA$year$b),nTST )
    
    
 ## Add location effect
   locHat=rep(NA, nTST )
   for(i in 1:nTST){   
            tmp<-grep(x=names(fm$ETA$loc$b),pattern=Y2$LOC[i])
            tmp<-ifelse(length(tmp)==0,0,fm$ETA$loc$b[tmp])
   			locHat[i]=tmp 
   	}

    
 ## Add genotype effect
   # Compute G21
   tmp<-which(colnames(X2)%in%colnames(X))
   X2=X2[,tmp]
   X2=scale(X2,scale=sdX,center=meansX)
   for(i in 1:ncol(X2)){
   	 tmp=is.na(X2[,i])
   	 if(any(tmp)){ X2[tmp,i]=0 }
   }
   G21<-tcrossprod(X2,X)/ncol(X)
   uHat=G21%*%GInv%*%PC.G%*%fm$ETA$G$b    
   gHat<-rep(NA,nTST) 
   for(i in 1:nTST){   
       tmp<-grep(x=rownames(X2),pattern=Y2$VAR[i])
       tmp<-ifelse(length(tmp)==0,0,uHat[tmp])
   	   gHat[i]=tmp 
   	}

   
 ## Add env. covariate effects
   tmp<-which(colnames(W2)%in%colnames(W))
   W2=W2[,tmp]
   W2=scale(W2,center=attributes(W)[[1]],scale=attributes(W)[[2]])
   W2<-W2/sqrt(ncol(W))
   wHat=W2%*%fm$ETA$W$b
 
   
 ## Add GW effect
  WW21<-tcrossprod(W2,W)
  ID.TST=as.integer(factor(x=Y2$VAR,levels=rownames(G21),ordered=T))
  ID.TRN= as.integer(factor(x=Y$VAR,levels=colnames(G21),ordered=T))
  GW21<-WW21*G21[ID.TST,ID.TRN]
  TMP=tcrossprod(GW21,GWInv)
  gwHat<-PC.GW%*%fm$ETA$GW$b
  gwHat<-TMP%*%gwHat

  yHat<- muHat + locHat   + gHat +  wHat  +  gwHat
  save(muHat,locHat,gHat,wHat,gwHat,file='predictions.RData')
  
timeEnd<-proc.time()
  
 ```
 
 getCor<-function(y,yHat,INDEX){
 	counts=table(INDEX)
 	COR<-rep(NA,nrow(counts))
 	for(i in 1:nrow(counts)){
 		tmp<-INDEX==names(counts)[i] ;
 		COR[i]=cor(y[tmp],yHat[tmp]) ;
 	}
    OUT=cbind(counts,COR)
    return(OUT)
 }

 COR=getCor(Y2$rdt,yHat,Y2$LOC)
 COR=COR[!is.na(COR[,2]),]
 
 pdf('correlations.pdf')
    plot(Y2$rdt~yHat,col=as.integer(factor(Y2$LOC)),cex=.5,ylab='Observed Yield',xlab='Predicted Yield',main=paste0('Correlation=',round(cor(yHat,Y2$rdt),3),'.'));
      abline(v=mean(yHat),lty=2,cex=2,col='white'); abline(v=mean(yHat),lty=2)
      abline(h=mean(Y2$rdt),lty=2,cex=2,col='white'); abline(h=mean(Y2$rdt),lty=2)
      abline(a=0,b=1,lwd=4,col='white');abline(a=0,b=1,col=1)      
      
     plot(x=COR[,1],y=COR[,2],xlab='Number of Records',ylab='Within Location Correlation',main=paste('Across Loc. Correlation=',round(cor(yHat,Y2$rdt),2)))
     lines(y=rep(mean(COR[COR[,1]>50,2]),2),x=c(50,150))
 dev.off()
 
 
 