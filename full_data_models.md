## Full-data  analysis

```R
 ## Parameters
  inputFile='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/standardized_data.RData'
  inputETA='/Users/gustavodeloscampos/Dropbox/arvalis/PIPELINES_2014/input/ETA.RData'
  outputFolder='/Users/gustavodeloscampos/WORK/ARVALIS/outputsGitHub/full_data_models/'

  library(BGLR)
  nIter=1200; burnIn=200
 ###
 
 dir.create(outputFolder) 
 setwd(outputFolder)
 load(inputFile)
 
 OUT=matrix(nrow=4,ncol=5,NA)
 colnames(OUT)=c('E','G','W','GxW','Error')
 rownames(OUT)=c('EL','EG','EGW','EGW_GxW')
 
 ### Model 1: Additive model (without markers)
  load(inputETA)
  ETA <- ETA$ETA1
  fmEL <- BGLR(y=Y$rdt,ETA=ETA,saveAt='EL_',nIter=nIter,burnIn=burnIn)
  OUT['EL','E']=fmEL$ETA$ENV$varB
  OUT['EL','G']=fmEL$ETA$VAR$varB
  OUT['EL','Error']=fmEL$varE
 
 ### Model 2: Additive model (with markers)
  load(inputETA)
  ETA <- ETA$ETA2
  fmEG <- BGLR(y=Y$rdt,ETA=ETA,saveAt='EG_',nIter=nIter,burnIn=burnIn)
  OUT['EG','E']=fmEG$ETA$ENV$varB
  OUT['EG','G']=fmEG$ETA$VAR$varB
  OUT['EG','Error']=fmEG$varE

 ### Model 3: Additive model (with markers and env. cov.)
  load(inputETA)
  ETA <- ETA$ETA3
  fmEGW <- BGLR(y=Y$rdt,ETA=ETA,saveAt='EGW_',nIter=nIter,burnIn=burnIn)
  OUT['EGW','E']=fmEGW$ETA$ENV$varB
  OUT['EGW','G']=fmEGW$ETA$VAR$varB
  OUT['EGW','W']=fmEGW$ETA$COV$varB
  OUT['EGW','Error']=fmEGW$varE
  
 ### Model 4: GxW Model 1 (Model 3 + interactions between markers and env. covariates)
  load(inputETA)
  ETA <- ETA$ETA4
  fmEGW_GxW <- BGLR(y=Y$rdt,ETA=ETA,saveAt='EGW_GxW_',nIter=nIter,burnIn=burnIn)
  OUT['EGW_GxW','E']=fmEGW_GxW$ETA$ENV$varB
  OUT['EGW_GxW','G']=fmEGW_GxW$ETA$VAR$varB
  OUT['EGW_GxW','W']=fmEGW_GxW$ETA$COV$varB
  OUT['EGW_GxW','GxW']=fmEGW_GxW$ETA$GxW$varB
  OUT['EGW_GxW','Error']=fmEGW_GxW$varE
  
  round(OUT,3)
```
 
Estimated Variace Components


| Model     | E    |  G     | W    | GxW   | Error |
| --------- |-----:| -----:| -----:| -----:| -----:| 
EL      | 201.068  | 14.644  |     --  |    --  | 22.638 | 
EG      |  200.821 |  16.937  |     --  |    --  | 22.659 | 
EGW      | 143.603 |  16.108 |  38.099 |     --  | 22.590 | 
EGW_GxW  | 144.390  | 14.044 |  34.216 |  5.990 |  18.311 | 


[Home](https://github.com/gdlc/ARVALIS/blob/master/README.md)
