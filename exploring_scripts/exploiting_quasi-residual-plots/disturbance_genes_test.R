# experiments noisy genes, to find sense for residual plots

# loading R function need ####
source("feature_code/R/VCR_pamr.R") #to import special vcr function made for pamr(NSC) classifier

# SRBCT dataset from the library pamr it's self ####
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
khanc$geneid=as.character(as.matrix(khan[,1]))
str(khanc) #63obs like in tibs 2002

#main body ####
#fitting
library(pamr)
?pamr.train

pamr=pamr.train(khanc) #data=SRBCT as pamr paper
new.scales=pamr.adaptthresh(pamr)
pamr=pamr.train(khanc, threshold.scale = new.scales)
pamr

?pamr.listgenes
pamr.listgenes(pamr, khanc, threshold = 3.7)

vcrpamr=vcr.pamr.train(data=khanc, pamrfit=pamr, threshold= 5)
vcrpamr$PAC

library(classmap)

qres=qresplot(vcrpamr$PAC, t(khanc$x)[,422], xlab = "Age (years)", opacity = 0.5,
         main = "quasi residual plot for male passengers",
         plotLoess = TRUE)


coeff=rep(NA,ncol(t(khanc$x)))

for (i in 1:length(coeff)) {
  lm=lm(vcrpamr$PAC~t(khanc$x)[,i])
  coeff[i]=lm$coefficients[2]
}

sorted_indices <- order(coeff, decreasing = TRUE)
sorted_indices[1:10]
which.max(coeff)
lm(vcrpamr$PAC~t(khanc$x)[,424])

#eliminating gene and seeing if confidence improve (this maybe can't be visible in analyzing only test error)
silplot(vcrpamr)
pamr.confusion(pamr, threshold=4.45)

classmap(vcrpamr, whichclass = 4)
