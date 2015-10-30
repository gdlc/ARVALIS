
BS<-function(x,knots=quantile(x,p=c(.2,.4,.6,.8)),degree=3,intercept=F){
	####
	# Construct basis functions for b-splines
	# should be equivalent to bs() of the spline package
	###
	n=length(x)
	X<-matrix(nrow=n,ncol=(degree-1+length(knots)),NA)
	for(i in 1:(degree-1)){
		X[,i]<-x^i
	}
	
	for(i in 1:length(knots)){
		cond<- x>knots[i]
		z<-(x-knots[i])^degree
		X[,i+degree-1]=z*cond
	
	}
	if(intercept){ X=cbind(1,X) }
	return(X)
	
}
