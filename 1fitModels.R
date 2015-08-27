rm(list=ls())
library(BGLR)
args=(commandArgs(TRUE))
if(length(args)==0){
	Mo="G_cholLW_GxW"
	thin=5
	colNum=0
  nIter<-30000
 burnIn<-10000
 }else{
	cat('eval agrs',args,"\n")
 for(i in 1:length(args)){
 	
    eval(parse(text=args[[i]]))
 }
}
cat("Mo:",Mo,"thin:",thin,"colNum:",colNum,"nIter:",nIter,"burnIn:",burnIn,"\n")
# Parameters
 dataDir<-file.path(Sys.getenv("USER_SCRATCH"),"ARVALIS/data")
 outDir<-file.path(Sys.getenv("USER_SCRATCH"),"ARVALIS/randReg/output",paste(nIter,"_",burnIn,sep=""))
 dir.create(outDir,recursive=T)
setwd(dataDir) 
load("DATA.rda")
library(splines)
library('BGLR')
 
 y<-Y$rdt
 LOCxYEAR=paste(Y$LOC,Y$YEAR,sep='-')
 VAR<-factor(x=Y$VAR,levels=rownames(X1),ordered=TRUE)
 Z.VAR<-as.matrix(model.matrix(~VAR-1))

 X<-cbind(X1,X2)
 G<-tcrossprod(X)/ncol(X)
 L<-t(chol(G))
 cholL.G=Z.VAR%*%L
 df0=.2 

setwd(outDir)
dir.create(Mo)
setwd(Mo)
#Models with all the EC and without the random regression on specific ec

#
if(Mo=="E_G_GxE"){
	load(file.path(dataDir,"cholL.GxE.rda"))
ETA=list(
        E=list(X=model.matrix(~factor(LOCxYEAR)-1),model='BRR'), 
		G=list(X=cholL.G,model="BRR"),
		GxE=list(X=cholL.GxE,model='BRR')
		)
 fm<-BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste(Mo,"_",sep=""))
 
 save(fm,file=paste('fm',Mo,'.rda',sep=""))
 }
 
#
if(Mo=="L_E_W"){
ETA=list(
		L=list(X=Z.VAR,model="BRR"),
		E=list(X=model.matrix(~factor(LOCxYEAR)-1),model='BRR'), 
		W=list(X=W/sqrt(ncol(W)),model='BRR')	
		)
 fm<-BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste(Mo,"_",sep=""))
 
 
 save(fm,file=paste('fm',Mo,'.rda',sep=""))
 }

# 
if(Mo=="G_W_GxW"){
	load(file.path(dataDir,"cholL.GxW.rda"))
ETA=list(
		G=list(X=cholL.G,model="BRR"),
		W=list(X=W/sqrt(ncol(W)),model='BRR'),
		GxW=list(X=cholL.GxW,model='BRR')
		)
 fm<-BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste(Mo,"_",sep=""))
 
 
 save(fm,file=paste('fm',Mo,'.rda',sep=""))
 }
#
if(Mo=="G_cholLW_GxW"){
	load(file.path(dataDir,"cholL.GxW.rda"))
	load(file.path(dataDir,"cholL.WW.rda"))
ETA=list(
        G=list(X=cholL.G,model="BRR"),
		W=list(X=cholL.WW,model='BRR'),
		GxW=list(X=cholL.GxW,model='BRR')
		)
 fm<-BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste(Mo,"_",sep=""))
 
 
 save(fm,file=paste('fm',Mo,'.rda',sep=""))
 } 
# 
if(Mo=="G_eigenW_GxW"){
	load(file.path(dataDir,"cholL.GxW.rda"))
	load(file.path(dataDir,"eigenWW.rda"))
ETA=list(
        G=list(X=cholL.G,model="BRR"),
		W=list(V=eigenWW$vectors,d=eigenWW$values,model='RKHS'),
		GxW=list(X=cholL.GxW,model='BRR')
		)
 fm<-BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste(Mo,"_",sep=""))
 
 
 save(fm,file=paste('fm',Mo,'.rda',sep=""))
 } 

#
if(Mo=="E_G_W_GxW"){
	load(file.path(dataDir,"cholL.GxW.rda"))
ETA=list(
        E=list(X=model.matrix(~factor(LOCxYEAR)-1),model='BRR'), 
		G=list(X=cholL.G,model="BRR"),
		W=list(X=W/sqrt(ncol(W)),model='BRR'),
		GxW=list(X=cholL.GxW,model='BRR')
		)
 fm<-BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste(Mo,"_",sep=""))
 
 
 save(fm,file=paste('fm',Mo,'.rda',sep=""))
 }
   
# Model with splines
if(Mo=="E_G_W_GxW_spF_spR" | Mo=="E_G_W_GxW_spF"){
 ec<-W[,colNum]
 ec.name=colnames(W)[colNum]
 dir.create(ec.name)
 setwd(ec.name)
 EC.ns=as.matrix(ns(ec,df=4))
 EC.ns<-scale(EC.ns,scale=TRUE,center=TRUE)
 x1=EC.ns[,1]
 x2=EC.ns[,2]
 x3=EC.ns[,3]
 x4=EC.ns[,4]

 # USING BGLR
 ECBS1cholL.G<-cholL.G; for(i in 1:nrow(cholL.G)){ ECBS1cholL.G[i,]<-x1[i]*cholL.G[i,] }
 ECBS2cholL.G<-cholL.G; for(i in 1:nrow(cholL.G)){ ECBS2cholL.G[i,]<-x2[i]*cholL.G[i,] }
 ECBS3cholL.G<-cholL.G; for(i in 1:nrow(cholL.G)){ ECBS3cholL.G[i,]<-x3[i]*cholL.G[i,] }
 ECBS4cholL.G<-cholL.G; for(i in 1:nrow(cholL.G)){ ECBS4cholL.G[i,]<-x4[i]*cholL.G[i,] } 
 W<-W[,-colNum]
 WW<-tcrossprod(W)/ncol(W)
GxW<-tcrossprod(cholL.G)*WW
diag(GxW)<-diag(GxW)+1/1000
cholL.GxW<-t(chol(GxW))

if(Mo=="E_G_W_GxW_spF"){
 ETA=list(
        E=list(X=model.matrix(~factor(LOCxYEAR)-1),model='BRR'), 
		G=list(X=cholL.G,model="BRR"),
		W=list(X=W/sqrt(ncol(W)),model='BRR'),
		GxW=list(X=cholL.GxW,model='BRR'),
		spF=list(X=EC.ns,model='FIXED')) 
}

if(Mo=="E_G_W_GxW_spF_spR"){
	ETA=list(
	    E=list(X=model.matrix(~factor(LOCxYEAR)-1),model='BRR'), 
		G=list(X=cholL.G,model="BRR"),
		W=list(X=W/sqrt(ncol(W)),model='BRR'),
		GxW=list(X=cholL.GxW,model='BRR'),
		spF=list(X=EC.ns,model='FIXED'),	
	    spR1=list(X=ECBS1cholL.G,model='BRR'),
		spR2=list(X=ECBS2cholL.G,model='BRR'),
		spR3=list(X=ECBS3cholL.G,model='BRR'),
		spR4=list(X=ECBS4cholL.G,model='BRR')
		)
    }
    
		 	 
 fm<-BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=paste(Mo,"_",sep=""))
 
 BHat2<-matrix(nrow=nrow(G),ncol=5,NA)
 BHat2[,1]<-L%*%fm$ETA$G$b+fm$mu
 BHat2[,2]<-fm$ETA$spF$b[1]+(if(Mo=="E_G_W_GxW_spF_spR") L%*%fm$ETA$spR1$b else 0)#do not use ifelse, because it will not return a vector, but only return the first element of vector
 BHat2[,3]<-fm$ETA$spF$b[2]+(if(Mo=="E_G_W_GxW_spF_spR") L%*%fm$ETA$spR2$b else 0)
 BHat2[,4]<-fm$ETA$spF$b[3]+(if(Mo=="E_G_W_GxW_spF_spR") L%*%fm$ETA$spR3$b else 0)
 BHat2[,5]<-fm$ETA$spF$b[4]+(if(Mo=="E_G_W_GxW_spF_spR") L%*%fm$ETA$spR4$b else 0)
 rownames(BHat2)<-rownames(G)

 YHat2<-matrix(nrow=length(ec),ncol=nrow(BHat2),0)
 for(i in 1:ncol(YHat2)){
   YHat2[,i]<-cbind(1,EC.ns)%*%BHat2[i,]
 }
 colnames(YHat2)=rownames(BHat2)
 counts<-table(Y$VAR)

 
 save(fm,YHat2,file='curves.RData')
 }
 