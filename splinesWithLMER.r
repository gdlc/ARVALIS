###
# Code for fitting, testing and ploting randome effects splines with lmer
###
#rm(list=ls())

## Change parameters here
 load('~/Dropbox/ARVALIS/data/Data2015/Y.RData')
 load('~/Dropbox/ARVALIS/data/Data2015/W.RData')
 minNEnv=20 # used to remove environments with few records.
 minNVAR=20 # used to remove genotypes with few records.
 nRecords=100 # used to determine which lines to include in the second plot.
 nPC=10 # number of env. covariates PC to be included.
 covName='Tmoyb0SEp1'
 colNum=which(colnames(W)==covName)
 LRT<-TRUE # perform Likelihood Ratio Test?
###

library(lme4)
library(splines)

for(i in 1:5){
  counts<-table(Y$VAR)
  tmp<-Y$VAR%in%names(counts)[which(counts>=minNVAR)]
  Y<-Y[tmp,]
  W<-W[tmp,]
  counts<-table(Y$ENV)
  tmp<-Y$ENV%in%names(counts)[which(counts>=minNEnv)] 
  Y<-Y[tmp,]
  W<-W[tmp,]
}


y<-Y$rdt
LOCxYEAR=paste(Y$LOC,Y$YEAR,sep='-')
VAR=Y$VAR
ec<-W[,colNum]
PC.W<-svd(scale(W[,-colNum]),nv=0,nu=nPC)$u

EC.ns=as.matrix(ns(ec,df=4,intercept=F))
EC.ns<-scale(EC.ns,scale=TRUE,center=TRUE)
x1=EC.ns[,1]
x2=EC.ns[,2]
x3=EC.ns[,3]
x4=EC.ns[,4]

fm0  <-lmer(y~PC.W+(1|VAR)+(1|LOCxYEAR),REML=F)
fm0    <-lmer(y~PC.W+(1|VAR)+(1|LOCxYEAR)+ec,REML=F)
fmFNS  <-lmer(y~PC.W+(1|VAR)+(1|LOCxYEAR)+EC.ns,REML=F)
fmRNS  <-lmer(y~PC.W+(1|LOCxYEAR)+x1+x2+x3+x4+(1|VAR)+(x1-1|VAR)+(x2-1|VAR)+(x3-1|VAR)+(x4-1|VAR),REML=F )

## Likelihood ratio test
if(LRT){
  AOV=matrix(nrow=3,ncol=5)
  rownames(AOV)<-c('Fixed-linear','Random-linear','Random-NS')
  colnames(AOV)<-c('DF','Deviance','FL','RL','RNS')
  comparison=anova(fm0,fmFNS,fmRNS)
  AOV[,1]=comparison$Df
  AOV[,2]=comparison$deviance
  AOV[1,4]= -(AOV[2,2]-AOV[1,2])
  AOV[1,5]= -(AOV[3,2]-AOV[1,2])
  AOV[2,5]= -(AOV[3,2]-AOV[2,2])
  AOV[2,3]=1-pchisq(AOV[1,2]-AOV[2,2],df=AOV[2,1]-AOV[1,1])
  AOV[3,3]=1-pchisq(AOV[1,2]-AOV[3,2],df=AOV[3,1]-AOV[1,1])
  AOV[3,4]=1-pchisq(AOV[2,2]-AOV[3,2],df=AOV[3,1]-AOV[2,1])
  print(AOV)
 }

BHat<-as.matrix(coef(fmRNS)$VAR)[,-c(2:(nPC+1))]
YHat<-matrix(nrow=length(ec),ncol=nrow(BHat),0)
for(i in 1:ncol(YHat)){
    YHat[,i]<-cbind(1,EC.ns)%*%BHat[i,]
}

yHat0<-cbind(1,EC.ns)%*%colMeans(BHat)

 
#pdf(paste0(covName,'.pdf'))
 plot(numeric()~numeric(),xlim=range(ec),ylim=range(YHat),xlab=covName,ylab='Predicted Yield',main='All')
 abline(v=ec,col=8,lty=1,lwd=.01)
 for(i in 1:ncol(YHat)){
    tmp<-ec[which(YHat[,i]==max(YHat[,i]))[1]]
    lines(x=sort(ec),y=YHat[order(ec),i],lwd=.3,col=2)
 }
 
 lines(x=sort(ec),y=yHat0[order(ec)],lwd=2,col=4)
 counts<-table(Y$VAR)
 
 IDs<-names(counts)[which(counts>nRecords)]
 plot(numeric()~numeric(),xlim=range(ec),ylim=range(YHat),xlab=covName,ylab='Predicted Yield',
      main=paste0('Linews with more than ',nRecords,'.'))
 abline(v=ec,col=8,lty=1,lwd=.01)
 for(i in 1:length(IDs)){
    newLine<-which(rownames(BHat)==IDs[i])
    lines(x=sort(ec),y=YHat[order(ec),newLine],col=2,lwd=.5)   
 }
 lines(x=sort(ec),y=yHat0[order(ec)],lwd=2,col=4)

#dev.off()


