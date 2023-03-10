# set session to project directory directory ####
library(here)
setwd(here())

# loading R function need ####
source("feature_code/R/VCR_pamr.R") #to import special vcr function made for pamr(NSC) classifier

# Main code of examples ####
## a generic microarray dataset form local ####
library(foreign) #to read WEKA arff format
datar=read.arff("sample_datasets/Ovarian.arff")
str(datar) #need manipulation to get to work with pamR
data=list() #pamr wants list
data$y=datar[,ncol(datar)]
data$x=datar[,1:(ncol(datar)-1)]
data$x=t(data$x) #pamr wants nvariable x obervations
str(data$x)
str(data$y)

#fitting
library(pamr)
?pamr.train

pamr=pamr.train(data) #data=SRBCT as pamr paper
pamr #chose a threshold by eye and get index
pamr$prob
yhat=pamr$yhat[3] #choose threshold 6.763
str(pamr$prob) #it's a 3d dim array [i,j,k]=[nobvs,class,thresholdindex]
pprob=pamr$prob[,,3]
ytrue=SRBCT$y

#producing the output for classmap
vcrpamr=vcr.pamr.train(data=data, pamrfit=pamr, threshold_index = 22) #data is feeded in same format that pam accepts

#silhouette visual plot
library(classmap)
?silplot #takes in a vcr out
silplot(vcrpamr) #classLabels = c("EWS","BL","NB","RMS") if you have names
pamr$nonzero[19]
pamr$se.scale
pamr$threshold.scale
pamr$call
pamr$threshold
pamr$scale.sd
str(pamr$centroids)
pamr$centroid.overall

#farness plot
classmap(vcrpamr, 2) #very very strange behaviour of the curve (opposite)



## on SRBCT dataset (Tibs 2002) ####
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
vcrpamr=vcr.pamr.train(data=SRBCT, pamrfit=pamr, threshold= 8) #data is feeded in same format that pam accepts
pamr.confusion(pamr, threshold=8)

#silhouette visual plot
library(classmap)
?silplot #takes in a vcr out
silplot(vcrpamr, classLabels = c("EWS","BL","NB","RMS") ) #classLabels = c("EWS","BL","NB","RMS") for SRBCT
pamr.confusion(pamr, threshold=8)  # check if there is match with class label (about line of ordering label)
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
classmap(vcrpamr, 3) #very very strange behaviour of the curve (opposite)

#### what can do
# check with other dataset
# try another farness (initfig) measures that you invent
# check all the process (all can be ok bu process can be flawed)
# try classic dataset

#troubleshooting
delta=(pamr$centroids-pamr$centroid.overall)/pamr$sd
delta.shrunk=soft.shrink(delta, pamr$threshold[19])
