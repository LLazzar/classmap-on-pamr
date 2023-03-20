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

# confusion ####
?pamr.confusion()
pamr.confusion(pamrtrain, threshold = 7.07)
pamrtrain$sample.subset



# cross-validation ####
?pamr.cv()
cvpamr=pamr.cv(pamrtrain, khanc)
cvpamr$folds #here we have indexes of each obvs left out of train for each set
             #crucial to reconstruct dd for pamr.cv
cvpamr$cv.objects[[1]] #if want to add pamr.cv as input of vcr.pamr shoulkd understand this object
cvpamr
# geneplot ####
pamr.geneplot(pamrtrain,khanc,threshold=7)
# indeterminate ####
pamr.indeterminate(pamrtrain$prob, mingap=0.2) #should use pamr predict
# listgenes ####
pamr.listgenes(pamrtrain, khanc, 4.2) #can be called only on pamr.train
# makeclass ####
pamr.makeclasses(khanc)
# menu ####


# plots ####
par(mar = c(6, 6, 6, 6) + 0.1)
pamr.plotcen(pamrtrain, khanc, 6)
pamr.plotcv(cvpamr)
dev.off()
pamr.plotcvprob(cvpamr,khanc,4) #the visualization we are trying to improve
pamr.plotfdr()
#
# predict ####
pamr.predict(pamrfit)


# pamr with survival data ####
gendata<-function(n=100, p=2000){
  tim <- 3*abs(rnorm(n))
  u<-runif(n,min(tim),max(tim))
  y<-pmin(tim,u)
  ic<-1*(tim<u)
  m <- median(tim)
  x<-matrix(rnorm(p*n),ncol=n)
  x[1:100, tim>m] <- x[1:100, tim>m]+3
  return(list(x=x,y=y,ic=ic))
}
# generate training data; 2000 genes, 100 samples
junk<-gendata(n=100)
y<-junk$y
ic<-junk$ic
x<-junk$x
d <- list(x=x,survival.time=y, censoring.status=ic,
          geneid=as.character(1:nrow(x)), genenames=paste("g",as.character(1:nrow(x)),sep=
                                                            ""))
# train model
a3<- pamr.train(d, ngroup.survival=4)
a3
pamr.plotstrata(a3, d$survival.time, d$censoring.status)
pamr.confusion.survival(a3, d$survival.time, d$censoring.status, a3$yhat[[3]])

