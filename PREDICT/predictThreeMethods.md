## This code illustrates predictions in testing sets using three different approaches

# (1) adding NAs in the entries for testing (fm0)
# (2) fitting the model to trn data and then deriving BLUP using standard regression (BLUP) formulae
# (3) same as 2 but using a reduced number of eigenvectors when fitting the model

library(BGLR)
data(wheat)
X=scale(wheat.X)
tst=1:100

y<-wheat.Y[,1]

G=tcrossprod(X)/ncol(X)

G11=G[-tst,-tst]
GInv=chol2inv(chol(G11))
G21=G[tst,-tst]


yNA=y
yNA[tst]=NA
fm0=BGLR(y=yNA,ETA=list(list(K=G,model='RKHS')),nIter=32000,burnIn=2000)
yHat0=fm0$yHat[tst]

fm=BGLR(y=y[-tst],ETA=list(list(K=G11,model='RKHS')),nIter=32000,burnIn=2000)
yHat<-G21%*%GInv%*%fm$ETA[[1]]$u

EVD=eigen(G11)

V=EVD$vectors[,EVD$values>.1]
d=EVD$values[EVD$values>.1]
fm2=BGLR(y=y[-tst],ETA=list(list(V=V,d=d,model='RKHS')),nIter=32000,burnIn=2000)
GInv2<-V%*%diag(1/d)%*%t(V)

yHat2<-G21%*%GInv%*%fm2$ETA[[1]]$u

cor(yHat0,yHat2)


unlink('*.dat')

plot(yHat0,yHat)



