## Cross validation 1 and 2

**Cross validation 1**.
Newly developed lines. Random assigns lines to folds. All the records for a given line are assigned to the same fold.
 
**Cross validation 2**.
To assess the ability of models to predict the performance of lines using the data collected in other environments.
 
```R
 ## Parameters
  inputFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/standardized_data.RData'
  inputETA='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/ETA.RData'
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/cross_validation/'
  
  #############
  #   CV    = 1:CV1,   2:CV2
  #   model = 1:EL,   2:EG,   3:EGW,   4:EGW_GxW
  #   seed  = 1 to 100
  #############
  CV <- 2
  model <- 3       
  nIter <- 20000
  burnIn <- 5000
  seed <- 1
 ###
 library(BGLR)

 dir.create(outputFolder) 
 setwd(outputFolder)
 load(inputFile)
 
 nfolds <- 5
 seeds <- sample(seq(1E4,1E5),100,replace=FALSE)
 models <- c("EL","EG","EGW","EGW_GxW")
 model <- models[model]
 y <- Y$rdt
 IDs <- rownames(X)
 IDy <- as.character(Y$VAR)
 folds <- rep(1:nfolds,each=ceiling(length(IDs)/nfolds))[order(runif(length(IDs)))]
 CVfolds <- rep(NA,nrow(Y))

 set.seed(seeds[seed])
 if(CV==1)
 {
   for(i in 1:nfolds)
   {   tmp <- IDy%in%IDs[folds==i]
       CVfolds[tmp] <- i
   }
 }

 if(CV==2)
 {
   for(i in IDs)
   {   tmp <- which(IDy==i)
       ni <- length(tmp)
       fold0 <- sample(1:nfolds,size=ni,replace=ni>nfolds)
       CVfolds[tmp] <- fold0
   }
 }

 ## Fitting model
 load(inputETA)
 ETA <- switch(model,
	'EL'      = ETA$ETA1, 
	'EG'      = ETA$ETA2,
	'EGW'     = ETA$ETA3,
	'EGW_GxW' = ETA$ETA4 		
 )
 yHatCV <- rep(NA,length(y))
 
 for(fold in 1:nfolds)
 {
   yNA <- y
   indexNA <- CVfolds==fold 
   yNA[indexNA] <- NA
   
   fm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn)
   yHatCV[indexNA] <- fm$yHat[indexNA]
 }

 OUT <- data.frame(y,yHatCV,ENV=Y$ENV)

 tt <- lapply(split(OUT,OUT$ENV),function(x)c(nrow(x),cor(x$y,x$yHatCV)))
 OUT2 <- do.call('rbind',tt)
 colnames(OUT2) <- c("n","corCV")
 OUT2 <- data.frame(ENV=rownames(OUT2),OUT2)

 ### Saving outputs 
 write.csv(OUT,paste0("Predictions_",model,"_CV",CV,".csv"),row.names=F)
 write.csv(OUT2,paste0("Predictions_by_Env_",model,"_CV",CV,".csv"),row.names=F)
```
