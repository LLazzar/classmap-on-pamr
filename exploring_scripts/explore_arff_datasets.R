library(foreign) #to read WEKA arff format
datar=read.arff("sample_datasets/arff_curda_rnaseq/GSE50760.arff")
str(datar) #need manipulation to get to work with pamR
dim(datar)
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

pamr.confusion(pamr,threshold = 10)


source("feature_code/R/VCR_pamr.R")

library(classmap)
vcrpamr=vcr.pamr.train(data,pamr, pamrfitcv=pamrcv, threshold=10)
vcrpamr$initfig


silplot(vcrpamr)
classmap(vcrpamr, whichclass =1, classLabels=c("1","2","3","4"), identify=TRUE, cutoff=0.5)


vcrpamrknn=vcr.knn.train(X=as.matrix(t(data$x)), as.factor(data$y), k=10)
silplot(vcrpamrknn)
classmap(vcrpamr, whichclass = 1)

library(randomForest)
forestfit=randomForest(x=as.matrix(t(data$x)), y=as.factor(data$y))
vcrpamrforest=vcr.forest.train(X=as.matrix(t(data$x)), y=as.factor(data$y), trainfit=forestfit)
silplot(vcrpamrforest)
?classmap()
