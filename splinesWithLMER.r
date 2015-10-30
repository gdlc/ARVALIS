###
# Code for fitting, testing and ploting randome effects splines with lmer
###

rm(list=ls())

## Change parameters here
 load('~/Dropbox/ARVALIS/data/Data2015/Y.RData')
 load('~/Dropbox/ARVALIS/data/Data2015/W.RData')
 nRecords=100 # used to determine which lines to plot
 covName='Tmoyb0SEp1'
 colNum=which(colnames(W)==covName)
###


library(lme4)
library(splines)

y<-Y$rdt
LOCxYEAR=paste(Y$LOC,Y$YEAR,sep='-')
VAR=Y$VAR
ec<-W[,colNum]
EC.ns=as.matrix(ns(ec,df=4))
EC.ns<-scale(EC.ns,scale=TRUE,center=TRUE)

x1=EC.ns[,1]
x2=EC.ns[,2]
x3=EC.ns[,3]
x4=EC.ns[,4]
fm0<-lmer(y~(1|VAR)+(1|LOCxYEAR))
 
fm0<-lmer(y~(1|VAR)+(1|LOCxYEAR))
fmL<-lmer(y~(1|VAR)+(1|LOCxYEAR)+ec)
fmNS<-lmer(y~(1|VAR)+(1|LOCxYEAR)+EC.ns)
fmRL<-lmer(y~(1|LOCxYEAR)+(1|VAR)+ec+(ec-1|VAR) )
fmRNS<-lmer(y~(1|LOCxYEAR)+x1+x2+x3+x4+(1|VAR)+(x1-1|VAR)+(x2-1|VAR)+(x3-1|VAR)+(x4-1|VAR) )
anova(fmRL,fmRNS)

BHat<-as.matrix(coef(fmRNS)$VAR)
YHat<-matrix(nrow=length(ec),ncol=nrow(BHat),0)
for(i in 1:ncol(YHat)){
    YHat[,i]<-cbind(1,EC.ns)%*%BHat[i,]
}

yHat0<-cbind(1,EC.ns)%*%colMeans(BHat)

 
pdf(paste0(covName,'.pdf'))
  plot(numeric()~numeric(),xlim=range(ec),ylim=range(YHat),xlab=covName,ylab='Predicted Yield',main='All')
  abline(v=ec,col=8,lty=1,lwd=.01)
  for(i in 1:ncol(YHat)){
     tmp<-ec[which(YHat[,i]==max(YHat[,i]))[1]]
     lines(x=sort(ec),y=YHat[order(ec),i],lwd=.3,col=2)
 }
 lines(x=sort(ec),y=yHat0[order(ec),lwd=2,col=4)
 counts<-table(Y$VAR)
 
 IDs<-names(counts)[which(counts>nRecords)]
 
 plot(numeric()~numeric(),xlim=range(ec),ylim=range(YHat),xlab=covName,ylab='Predicted Yield',main=paste0('Linews with more than ',nRecords,'.'))
 abline(v=ec,col=8,lty=1,lwd=.01)
 
 for(i in 1:length(IDs)){
    newLine<-which(rownames(BHat)==IDs[i])
    lines(x=sort(ec),y=YHat[order(ec),newLine],col=2,lwd=.5)   
 }
 lines(x=sort(ec),y=yHat0[order(ec)],lwd=2,col=4)

dev.off()