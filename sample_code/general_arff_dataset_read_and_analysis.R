library(foreign) #to read WEKA arff format
datar=read.arff("sample_datasets/Breast_GSE26304.arff")
str(datar) #need manipulation to get to work with pamR
data=list() #pamr wants list
data$y=datar[,ncol(datar)]
data$x=datar[,1:(ncol(datar)-1)]
data$x=t(data$x) #pamr wants nvariable x obervations
str(data$x)
str(data$y)

library(pamr)

pamr=pamr.train(data)
threshold.scale=pamr.adaptthresh(pamr)

pamr=pamr.train(data, threshold.scale = threshold.scale)
pamr

pamrcv=pamr.cv(pamr, data) #best thres 1.6
pamrcv

pamr.confusion(pamr,threshold = 1.2)


source("feature_code/R/VCR_pamr.R")

library(classmap)
vcrpamr=vcr.pamr.train(data,pamr, pamrfitcv=pamrcv, threshold=1.6)

silplot(vcrpamr)
classmap(vcrpamr, whichclass = 4)
