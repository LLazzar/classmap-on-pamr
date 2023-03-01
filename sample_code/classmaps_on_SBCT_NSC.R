#set session to current file dir
#setwd("C:/Users/Admin/OneDrive - Alma Mater Studiorum Università di Bologna/maybethesis/classmap_on_pamr")

source("VCR_auxiliaryFunctions.R") #to import functions needed for VCR_pamr (especially checkLabels)
source("VCR_pamr.R") #to import special vcr function made for pamr(NSC) classifier
library(cellWise) #for transfo function used in Comp fareness in VCR-auxillary

#load microarray generic dataset from local file
library(foreign) #to read WEKA arff format
datar=read.arff("Ovarian.arff")
str(datar) #need manipulation to get to work with pamR
data=list()
data$y=datar[,ncol(datar)]
data$x=datar[,1:(ncol(datar)-1)]
data$x=t(data$x)
str(data$x)
str(data$y)

#load SBCT dataset from library plsgenomics
library(plsgenomics)
data(SRBCT)
?SRBCT
str(SRBCT) #already ready for pamr, only renames X->x, Y->y
SRBCT$x=t(SRBCT$X) #pamr wants genes(rows)xsamples(columns)
SRBCT$y=SRBCT$Y

#fitting
library(pamr)
?pamr.train

pamr=pamr.train(SRBCT) #data=SRBCT as pamr paper
pamr #threshold 6.763 seems better, correspond to index 19
yhat=pamr$yhat[3] #choose threshold 6.763
str(pamr$prob) #it's a 3d dim array [i,j,k]=[nobvs,class,thresholdindex]
pprob=pamr$prob[,,3]
ytrue=SRBCT$y

#producing the output
vcrpamr=vcr.pamr.train(data=SRBCT, pamrfit=pamr, threshold_index = 19) #data is feeded in same format that pam accepts

#silhouette visual plot
library(classmap)
?silplot #takes in a vcr out
silplot(vcrpamr, classLabels = c("EWS","BL","NB","RMS") ) #classLabels = c("EWS","BL","NB","RMS") for SRBCT
pamr$nonzero[19]
pamr$barbara
pamr$se.scale
pamr$threshold.scale
pamr$call
pamr$threshold
pamr$scale.sd
str(pamr$centroids)
pamr$centroid.overall

#farness plot
classmap(vcrpamr, 4) #very very strange behaviour of the curve (opposite)

#### what can do 
# check with other dataset
# try another farness (initfig) measures that you invent
# check all the process (all can be ok bu process can be flawed)
# try classic dataset

#troubleshooting 
delta=(pamr$centroids-pamr$centroid.overall)/pamr$sd
delta.shrunk=soft.shrink(delta, pamr$threshold[19])
