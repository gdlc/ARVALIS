###
# Code for fitting, testing and ploting random effects splines with lmer
##

rm(list=ls())
 library(lme4)
 library(splines)
 
## Change parameters here
 phenotypeFile<-'~/Dropbox/ARVALIS/data/Data2015/Y.RData'
 envCovFile<-'~/Dropbox/ARVALIS/data/Data2015/W.RData'
 nPC=10 # number of env. covariates PC to be included.
 outputFolder='~/Dropbox/ARVALIS/SPLINES/output'
 covNum<-NULL # by default it runs all covariates, otherwise specify a vector in integers
 plotLines<-FALSE
 nRecordsPrint<-100
###
 
 
 dir.create(outputFolder)
 setwd(outputFolder)
 load(phenotypeFile)
 load(envCovFile)
 
 tmp<-names(table(Y$VAR))[table(Y$VAR)>100]
 W<-W[Y$VAR%in%tmp,]
 Y<-Y[Y$VAR%in%tmp,]


 y<-Y$rdt
 LOCxYEAR=paste(Y$LOC,Y$YEAR,sep='-')
 VAR=Y$VAR
 
  counts<-table(Y$VAR)
 IDs<-names(counts)[which(counts>nRecordsPrint)]

 covNames<-colnames(W)
 pdf('plots.pdf')
 fileOut<-file('logLik_varE.txt',open='w') 
 msg0<-c( c('ll_Int','ll_LiF' ,'ll_SpF' ,'ll_LiR' ,'ll_SpR'),
 		 c('BIC_LiF','BIC_SpF','BIC_LiR','BIC_SpR'),
 		 'varExp'
 	)
 write(file=fileOut,append=T,x=msg0,ncol=length(msg0))
 if(is.null(covNum)){ covNum=1:ncol(W) } 
 W=scale(W)
 

 for(i in covNum){

  covName<-covNames[i]
  cat('Working cov #',i,'(',covName,')','\n')
  ec<-W[,i]
  ec<-ec+runif(min=-sd(ec)/1000,max=sd(ec)/1000,n=length(ec))
  EC.ns=ns(ec,df=4)
  
  varX=apply(MARGIN=2,X=EC.ns,FUN=var)  
  x1=EC.ns[,1]
  x2=EC.ns[,2]
  x3=EC.ns[,3]
  x4=EC.ns[,4]

  PC.W<-svd(W[,-i],nv=0,nu=nPC)$u
  PC.W<-scale(PC.W,scale=F,center=T)
  fmInt <-lmer(y~(1|VAR)+(1|LOCxYEAR)+PC.W,REML=F)
 
  fmLiF <-lmer(y~(1|VAR)+(1|LOCxYEAR)+PC.W+ec,REML=F)
  fmSpF <-lmer(y~(1|VAR)+(1|LOCxYEAR)+PC.W+x1+x2+x3+x4,REML=F) 

  fmLiR <-lmer(y~(1|VAR)+(1|LOCxYEAR)+PC.W+ec+(ec-1|VAR),REML=F)
  fmSpR<-lmer( y~(1|VAR)+(1|LOCxYEAR)+PC.W+x1+x2+x3+x4+(x1-1|VAR)+(x2-1|VAR)+(x3-1|VAR)+(x4-1|VAR),REML=F )
  
  varB<-as.data.frame(VarCorr(fmSpR))[as.data.frame(VarCorr(fmSpR))[,2]%in%c('x1','x2','x3','x4'),4][4:1]
  varExp<-crossprod(varX,varB)
  BICInt<-BIC(fmInt)
  msg<-c( logLik(fmInt),logLik(fmLiF),logLik(fmSpF),logLik(fmLiR),logLik(fmSpR),
  	      BICInt-BIC(fmLiF),BICInt-BIC(fmSpF),BICInt-BIC(fmLiR),BICInt-BIC(fmSpR),
  	      varExp
  	    )
  names(msg)=msg0
  print(round(msg,2))      
  
  write(file=fileOut,append=T,x=msg,ncol=length(msg))

  BHat=as.matrix((coef(fmSpR)$VAR[,c('(Intercept)','x1','x2','x3','x4')]))
  YHat<-matrix(nrow=length(ec),ncol=nrow(BHat),0)
  
  for(j in 1:ncol(YHat)){
    YHat[,j]<-cbind(1,EC.ns)%*%BHat[j,]
  }
  
  plot(numeric()~numeric(),xlim=range(ec),ylim=range(YHat),xlab=covName,ylab='Predicted Yield',
      main=paste0('Lines with more than ',nRecordsPrint,' records.'))
 	abline(v=ec,col=8,lty=1,lwd=.01)
 	for(i in 1:length(IDs)){
    	newLine<-which(rownames(BHat)==IDs[i])
    	lines(x=sort(ec),y=YHat[order(ec),newLine],col=2,lwd=.5)   
 	}
 	lines(x=sort(ec),y=rowMeans(YHat)[order(ec)],col=4,lwd=2)
}
dev.off()
close(fileOut)







  # Plots 
  
    BHatInt<-as.matrix(coef(fmInt)$VAR)[,1]
  BHatLinear<-as.matrix(coef(fmLinear)$VAR)[,2:1]
  BHatSpline<-as.matrix(coef(fmSpline)$VAR)[,1:5]

  YHatInt<-matrix(nrow=length(ec),ncol=nrow(BHatLinear),0)
  YHatLinear<-YHatInt
  YHatSpline<-YHatInt
  
  
  for(i in 1:ncol(YHatInt)){
    YHatSpline[,i]<-cbind(1,EC.ns)%*%BHatSpline[i,]
  }
  

  if(plotLines){

		plot(numeric()~numeric(),xlim=range(ec),ylim=range(YHat),xlab=covName,ylab='Predicted Yield',
      main=paste0('Lines with more than ',nRecordsPrint,' records.'))
 	abline(v=ec,col=8,lty=1,lwd=.01)
 	for(i in 1:length(IDs)){
    	newLine<-which(rownames(BHat)==IDs[i])
    	lines(x=sort(ec),y=YHat[order(ec),newLine],col=2,lwd=.5)   
 	}
 	lines(x=sort(ec),y=rowMeans(YHat)[order(ec)],col=4,lwd=2)
 	
   }
   # end plot NS
   



TMP=read.table('logLik_varE.txt',header=T,stringsAsFactors=F)


## 
 significance=.01
 nTests=nrow(TMP)
 threshold=-log10(significance/nTests)
 
 
## Has the covariate a linear association?
 lrtStat=-2*(TMP[,'ll_Int']-TMP[,'ll_LiF'])
 DF=1
 negLogPValue=-log10(pchisq(q=lrtStat,df=DF,lower.tail=F))
 plot(negLogPValue,cex=.4,col=8, main='H0: no covariate, Ha: linear-fixed',ylab='-log10(p-value)')
  abline(h=threshold)
  tmp<-which(negLogPValue>threshold)
  text(labels=colnames(W)[tmp],x=tmp,y=negLogPValue[tmp],cex=.6)

  
## Non-linear association?
 lrtStat=-2*(TMP[,'ll_Int']-TMP[,'ll_SpF'])
 DF=4
 negLogPValue=-log10(pchisq(q=lrtStat,df=DF,lower.tail=F))
 plot(negLogPValue,cex=.2,col=8, main='H0: no covariate, Ha: spline-fixed',ylab='-log10(p-value)')
  abline(h=threshold)
  tmp<-which(negLogPValue>threshold)
  text(labels=colnames(W)[tmp],x=tmp,y=negLogPValue[tmp],cex=.6)  
  
  
## Beyond linear?
 lrtStat=-2*(TMP[,2]-TMP[,3])
 DF=3
 negLogPValue=-log10(pchisq(q=lrtStat,df=DF,lower.tail=F))
 plot(negLogPValue,cex=.2,col=8, main='H0: linear-fixed, Ha: spline-fixed',ylab='-log10(p-value)')
  abline(h=threshold)
  tmp<-which(negLogPValue>threshold)
  text(labels=colnames(W)[tmp],x=tmp,y=negLogPValue[tmp],cex=.6)  
  
## GxE? (linear)
	lrtStat=-2*(TMP[,'ll_LiF']-TMP[,'ll_LiR'])
	DF=1
	negLogPValue=-log10(pchisq(q=lrtStat,df=DF,lower.tail=F))
 	plot(negLogPValue,cex=.2,col=8, main='H0: linear-fixed, Ha: linear-random',ylab='-log10(p-value)')
  	abline(h=threshold)
  	tmp<-which(negLogPValue>threshold)
  	text(labels=colnames(W)[tmp],x=tmp,y=negLogPValue[tmp],cex=.6)  

## GxE? (SPLINE)
	lrtStat=-2*(TMP[,'ll_SpF']-TMP[,'ll_SpR'])
	DF=4
	negLogPValue=-log10(pchisq(q=lrtStat,df=DF,lower.tail=F))
 	plot(negLogPValue,cex=.2,col=8, main='H0: spline-fixed, Ha: spline-random',ylab='-log10(p-value)')
  	abline(h=threshold)
  	tmp<-which(negLogPValue>threshold)	
	text(labels=colnames(W)[tmp],x=tmp,y=negLogPValue[tmp],cex=.6)  


## Error Variance
	
	propVar=(TMP[,'var_Int']-TMP[,'var_LiF'])/TMP[,'var_Int']

	propVar=(TMP[,'var_Int']-TMP[,'var_SpF'])/TMP[,'var_Int']

	propVar=(TMP[,'var_Int']-TMP[,'var_LiR'])/TMP[,'var_Int']

	propVar=(TMP[,'var_Int']-TMP[,'var_SpR'])/TMP[,'var_Int']

#######################################
#  Code below here is not being used  #
#######################################

if(FALSE){

 if(plotLines){ pdf('plotd.pdf') }

  -log10(anova(fmInt,fmLiF)[8][[1]][2]), -log10(anova(fmInt,fmSpF)[8][[1]][2]),-log10(anova(fmLiF,fmLiR)[8][[1]][2]),-log10(anova(fmSpF,fmSpR)[8][[1]][2]) )
 

 
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
  

   
  
  -log10(anova(fmInt,fmLinear)[8][[1]][2]),-log10(anova(fmInt,fmSpline)[8][[1]][2]),-log10(anova(fmLinear,fmSpline)[8][[1]][2]))
  msg<-c(msg,cor(as.vector(YHatInt),as.vector(YHatSpline))^2,cor(as.vector(YHatLinear),as.vector(YHatSpline))^2)
  names(msg)<-c('LRT_IntVsLinear','LRT_IntVsSpline','LRT_LinearVsSpline','R2_Int','R2_Linear')
  print(msg)
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
 
 
 	lines(x=sort(ec),y=yHat0[order(ec)],lwd=2,col=4)
   }
   # end plot NS
   
}
close(fileLRT_R2)

dev.off()


