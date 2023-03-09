#SRBCT dataset from library genomic
library(plsgenomics)
data(SRBCT)
?SRBCT
str(SRBCT) #already ready for pamr, only renames X->x, Y->y
SRBCT$x=t(SRBCT$X) #pamr wants genes(rows)xsamples(columns)
SRBCT$y=SRBCT$Y
str(SRBCT)

#SRBCT dataset from the library pamr it's self
data(khan) #used in tibs 2002 PNAS
str(khan) #still needs preprocessing
length(khan)
ncol(khan)
khanc=list()
khanc$x=khan[2:nrow(khan),3:ncol(khan)] #we have factors to be converted in numeric

for (col in names(khanc$x)[sapply(khanc$x, is.factor)]) {
  khanc$x[[col]] <- as.numeric(as.character(khanc$x[[col]]))
}

khanc$x=as.matrix(khanc$x)
str(khanc$x) #now ok

khanc$y=as.factor(as.matrix(khan[1,3:ncol(khan)]))
str(khanc$y)
#khanc$geneid=khan[,2]
str(khanc) #63obs like in tibs 2002


#fitting and exploring
library(pamr)
?pamr.train

pamrtrain=pamr.train(khanc) #basic fit
pamrtrain

# testing result of applying classmap ####
source("feature_code/R/VCR_pamr.R")
library(classmap)
vcrpamr=vcr.pamr.train(data=khanc, pamrfit=pamrtrain, threshold_index = 17)
silplot(vcrpamr) #roughly like same info of figure 5 of tibs 2002

# adaptive thresholds ####
thresh=pamr.adaptthresh(pamrtrain) #producing thresholds scale different for class
pamrtrain.adapt=pamr.train(khanc, threshold.scale = thresh )
pamrtrain.adapt

# pamr.confusion ####
?pamr.confusion()
pamr.confusion(pamrtrain, threshold = 7.07)
pamrtrain$sample.subset



