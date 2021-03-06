## Full-data  analysis

```R
 ## Parameters
  inputFolder='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/'
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/full_data_models/'

  nIter=1200; burnIn=200 #use many more iterations, this is just for testing the code.
 ###
 library(BGLR)
 
 dir.create(outputFolder) 
 setwd(outputFolder)
 load(paste0(inputFolder,"/standardized_data.RData"))
 
 OUT=matrix(nrow=4,ncol=5,NA)
 colnames(OUT)=c('E','G','W','GxW','Error')
 rownames(OUT)=c('EL','EG','EGW','EGW_GxW')
 
 ### Model 1: Additive model (without markers)
  load(paste0(inputFolder,"/ETA_EL.RData"))
  fmEL <- BGLR(y=Y$y,ETA=ETA,saveAt='EL_',nIter=nIter,burnIn=burnIn)
  OUT['EL','E']=fmEL$ETA$ENV$varB
  OUT['EL','G']=fmEL$ETA$VAR$varB
  OUT['EL','Error']=fmEL$varE
 
 ### Model 2: Additive model (with markers)
  load(paste0(inputFolder,"/ETA_EG.RData"))
  fmEG <- BGLR(y=Y$y,ETA=ETA,saveAt='EG_',nIter=nIter,burnIn=burnIn)
  OUT['EG','E']=fmEG$ETA$ENV$varB
  OUT['EG','G']=fmEG$ETA$VAR$varB
  OUT['EG','Error']=fmEG$varE

 ### Model 3: Additive model (with markers and env. cov.)
  load(paste0(inputFolder,"/ETA_EGW.RData"))
  fmEGW <- BGLR(y=Y$y,ETA=ETA,saveAt='EGW_',nIter=nIter,burnIn=burnIn)
  OUT['EGW','E']=fmEGW$ETA$ENV$varB
  OUT['EGW','G']=fmEGW$ETA$VAR$varB
  OUT['EGW','W']=fmEGW$ETA$COV$varB
  OUT['EGW','Error']=fmEGW$varE
  
 ### Model 4: GxW Model 1 (Model 3 + interactions between markers and env. covariates)
  load(paste0(inputFolder,"/ETA_EGW_GxW.RData"))
  fmEGW_GxW <- BGLR(y=Y$y,ETA=ETA,saveAt='EGW_GxW_',nIter=nIter,burnIn=burnIn)
  OUT['EGW_GxW','E']=fmEGW_GxW$ETA$ENV$varB
  OUT['EGW_GxW','G']=fmEGW_GxW$ETA$VAR$varB
  OUT['EGW_GxW','W']=fmEGW_GxW$ETA$COV$varB
  OUT['EGW_GxW','GxW']=fmEGW_GxW$ETA$GxW$varB
  OUT['EGW_GxW','Error']=fmEGW_GxW$varE
  
  round(OUT,3)
```
 
Estimated Variance Components (Old Data)


| Model     | E    |  G     | W    | GxW   | Error |
| --------- |-----:| -----:| -----:| -----:| -----:| 
EL      | 201.068  | 14.644  |     --  |    --  | 22.638 | 
EG      |  200.821 |  16.937  |     --  |    --  | 22.659 | 
EGW      | 143.603 |  16.108 |  38.099 |     --  | 22.590 | 
EGW_GxW  | 144.390  | 14.044 |  34.216 |  5.990 |  18.311 | 


 
Estimated Variance Components (Data 2015)


| Model     | E    |  G     | W    | GxW   | Error |
| --------- |-----:| -----:| -----:| -----:| -----:| 
EL      | 201.068  | 14.644  |     --  |    --  | 22.638 | 
EG      |  200.821 |  16.937  |     --  |    --  | 22.659 | 
EGW      | 143.603 |  16.108 |  38.099 |     --  | 22.590 | 
EGW_GxW  | 144.390  | 14.044 |  34.216 |  5.990 |  18.311 | 

[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
