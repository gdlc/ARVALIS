###
# Code for fitting, testing and ploting randome effects splines with lmer
###
#  Questions
#        - Rnadom intercept: basic model heritability
#        - Random slope LRT randInt Vs randReg & R-sq
#        - Random shape LRT randReg Vs ranSpline
#        -
# 119 se cae y sigue cayendose 120, 121....
#
rm(list=ls())
 library(lme4)
 library(splines)

 getR2<-function(YHat,YHat0){
    tmp<-matrix(NA,nrow=ncol(YHat),ncol=ncol(YHat0))
    for(i in 1:ncol(YHat)){ 
       for(j in 1:ncol(YHat0)){
         tmp[i,j]<-cor(YHat[,i],YHat0[,j])^2
       }
       
    }
    out= colMeans(tmp)
    return(out)
    
 }
 
 
## Change parameters here
 phenotypeFile<-'~/Dropbox/ARVALIS/data/Data2015/Y.RData'
 envCovFile<-'~/Dropbox/ARVALIS/data/Data2015/W.RData'
 nRecordsPrint=100 # used to determine which lines to include in the second plot.
 nPC=10 # number of env. covariates PC to be included.
 outputFolder='~/tmp_arvalis/'
 covNum<-NULL # by default it runs all covariates, otherwise specify a vector in integers
 plotLines<-FALSE
###


 setwd(outputFolder)
 load(phenotypeFile)
 load(envCovFile)
 y<-Y$rdt
 LOCxYEAR=paste(Y$LOC,Y$YEAR,sep='-')
 VAR=Y$VAR
 
 if(plotLines){ pdf('plotd.pdf') }
 
 covNames<-colnames(W)
 fileLRT_R2<-file('LRT.txt',open='w') 
 msg<-c('RandInt_RandLinear','RandInt_RandSpline','RandLinear_RandSpline','R2_Linear','R2_Spline')
 write(file=fileLRT_R2,append=T,x=msg,ncol=length(msg))
 if(is.null(covNum)){ covNum=1:ncol(W) } 
 W=scale(W)
 for(i in 119:123){ W[,i]=W[,i]+runif(min=-.001,max=.001,n=nrow(W))}
 
 for(i in covNum){
  covName<-covNames[i]
  cat('Working cov #',i,'(',covName,')','\n')
  ec<-W[,i]
  if(nPC>0){PC.W<-svd(W[,-i],nv=0,nu=nPC)$u}
  
  # Natural spline
  EC.ns=as.matrix(ns(ec,df=4,intercept=F))
  EC.ns<-scale(EC.ns,scale=TRUE,center=TRUE)
  x1=EC.ns[,1]
  x2=EC.ns[,2]
  x3=EC.ns[,3]
  x4=EC.ns[,4]

  if(nPC>0){
    fmInt <-lmer(y~(1|VAR)+(1|LOCxYEAR)+x1+x2+x3+x4+PC.W,REML=F)
    fmLinear <-lmer(y~(1|VAR)+(1|LOCxYEAR)+(ec-1|VAR)+x1+x2+x3+x4+PC.W,REML=F)
    fmSpline<-lmer(y~(1|LOCxYEAR)+x1+x2+x3+x4+(1|VAR)+(x1-1|VAR)+(x2-1|VAR)+(x3-1|VAR)+(x4-1|VAR)+PC.W,REML=F )
  }else{
    fmInt <-lmer(y~(1|VAR)+(1|LOCxYEAR)+x1+x2+x3+x4,REML=F)
    fmLinear <-lmer(y~(1|VAR)+(1|LOCxYEAR)+(ec-1|VAR)+x1+x2+x3+x4,REML=F)
    fmSpline<-lmer(y~(1|LOCxYEAR)+x1+x2+x3+x4+(1|VAR)+(x1-1|VAR)+(x2-1|VAR)+(x3-1|VAR)+(x4-1|VAR),REML=F )
  }

 
  BHatInt<-as.matrix(coef(fmInt)$VAR)[,1]
  BHatLinear<-as.matrix(coef(fmLinear)$VAR)[,2:1]
  BHatSpline<-as.matrix(coef(fmSpline)$VAR)[,1:5]

  YHatInt<-matrix(nrow=length(ec),ncol=nrow(BHatLinear),0)
  YHatLinear<-YHatInt
  YHatSpline<-YHatInt
  
  
  for(i in 1:ncol(YHatInt)){
    YHatInt[,i]<-BHatInt[i]
    YHatLinear[,i]<-cbind(1,ec)%*%BHatLinear[i,]
    YHatSpline[,i]<-cbind(1,EC.ns)%*%BHatSpline[i,]
  }
  
  
  yHat0<-cbind(1,EC.ns)%*%colMeans(BHat)

  msg<-c( -log10(anova(fmL.F,fmL.R)[8][[1]][2]),-log10(anova(fmL.F,fmNS.R)[8][[1]][2]),-log10(anova(fmL.R,fmNS.R)[8][[1]][2]))
  msg<-c(msg,getR2(YHat,cbind(ec,yHat0)))
  msg=round(msg,4)
  
  write(file=fileLRT_R2,append=T,x=msg,ncol=length(msg))

  # Plot of random NS
  if(plotLines){
 	plot(numeric()~numeric(),xlim=range(ec),ylim=range(YHat),xlab=covName,ylab='Predicted Yield',main='All')
 	abline(v=ec,col=8,lty=1,lwd=.01)
 	for(i in 1:ncol(YHat)){
    	tmp<-ec[which(YHat[,i]==max(YHat[,i]))[1]]
    	lines(x=sort(ec),y=YHat[order(ec),i],lwd=.3,col=2)
 	}
 
 	lines(x=sort(ec),y=yHat0[order(ec)],lwd=2,col=4)
 	counts<-table(Y$VAR)
 
 	IDs<-names(counts)[which(counts>nRecordsPrint)]
 	plot(numeric()~numeric(),xlim=range(ec),ylim=range(YHat),xlab=covName,ylab='Predicted Yield',
      main=paste0('Lines with more than ',nRecordsPrint,' records.'))
 	abline(v=ec,col=8,lty=1,lwd=.01)
 	for(i in 1:length(IDs)){
    	newLine<-which(rownames(BHat)==IDs[i])
    	lines(x=sort(ec),y=YHat[order(ec),newLine],col=2,lwd=.5)   
 	}
 	lines(x=sort(ec),y=yHat0[order(ec)],lwd=2,col=4)
   }
   # end plot NS
   
}
close(fileLRT_R2)

dev.off()


