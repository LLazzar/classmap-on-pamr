# Preparing the gene-expr microarray dataset SBCT

library(foreign) #to read WEKA arff format
datar=read.arff("sample_datasets/SRBCT.arff")
str(datar) #need manipulation to get to work with pamR
train_indices <- sample(nrow(datar), nrow(datar) * 0.8)  # 70% for training
traindata=list() #pamr wants list
traindata$y=datar[train_indices,ncol(datar)]
traindata$x=datar[train_indices,1:(ncol(datar)-1)]
traindata$x=t(traindata$x) #pamr wants nvariable x obervations
str(traindata$x)
str(traindata$y)
testdata=list()
testdata$y=datar[-train_indices,ncol(datar)]
testdata$x=datar[-train_indices,1:(ncol(datar)-1)]
testdata$x=t(testdata$x) #pamr wants nvariable x obervations
str(testdata$x)
str(testdata$y)

# Fitting
library(pamr)
?pamr.train
pamr=pamr.train(traindata) #data=SRBCT as pamr paper
pamr #chose a threshold by eye and get index
pamr$prob
yhat=pamr$yhat[3] #choose threshold 6.763
str(pamr$prob) #it's a 3d dim array [i,j,k]=[nobvs,class,thresholdindex]
pprob=pamr$prob[,,3]
ytrue=SRBCT$y
pamrcv=pamr.cv(pamr,traindata)
pamrcv

# Visualizing on training set though classmap

source("feature_code/R/VCR_pamr.R") #to import special vcr function made for pamr(NSC) classifier
vcrpamr=vcr.pamr.train(data=traindata, pamrfit=pamr, threshold = 6) #data is feeded in same format that pam accepts
vcrpamrcv=vcr.pamr.train(data=traindata, pamrfit=pamr, pamrfitcv=pamrcv, threshold = 6)

## silhouette visual plot
library(classmap)
?silplot #takes in a vcr out
silplot(vcrpamr) #classLabels = c("EWS","BL","NB","RMS") if you have names
silplot(vcrpamrcv)

## classmap/farness plot visual
?classmap
classmap(vcrpamr, whichclass=4)


# Predicting test set
?pamr.predict
ypred=pamr.predict(fit=pamr, newx=testdata$x, threshold=6, type="class")
table(ypred,testdata$y) #all correct
testpost=pamr.predict(fit=pamr, newx=testdata$x, type="posterior") #not needed

#Visualizing on the test set though classmap
vcr.pamr.newdata(newdata = testdata, vcr.pamr.train.out = vcrpamr )



########################################################

labels=unlist(pamr$yhat[20], use.names = FALSE)
colors <- c("red", "green", "cyan","purple")
col_vec <- colors[match(labels, c("  1", "  2", "  3", "  4"))]
labels=as.numeric(as.factor(data$y))
col_vec_true = colors[match(data$y, c("  1", "  2", "  3", "  4"))]

intens=rep(NA,length(col_vec_true))
for (i in 1:length(col_vec_true)){
  pal <- colorRamp(c(col_vec_true[i], "white"))
  intens[i] = rgb(pal(0.9*vcrpamr$PAC[i]), maxColorValue = 255)
}

posid=vcrpamr$posid
sd=vcrpamr$sd
n=ncol(data$x)
sdr=sd[posid]
xr=data$x[posid,]
pairwiser=matrix(NA, nrow = n, ncol = n)
xr<- t(apply(t(xr), 1, function(x) x / sdr))
pairwiser=dist(xr)

sd=vcrpamr$sd
n=ncol(data$x)
pairwise=matrix(NA, nrow = n, ncol = n)
x<- t(apply(t(data$x), 1, function(x) x / sd))
pairwise=dist(x)

mds=cmdscale(pairwiser, k=2)
plot(mds, col=col_vec_true, pch=21, bg=intens, cex=0.8)


##########################################################


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

pamr=pamr.train(SRBCT ,gene.subset = c(1:2300)) #data=SRBCT as pamr paper
pamr
cvpamr=pamr.cv(pamr, SRBCT)
cvpamr
pamr #threshold 6.763 seems better, correspond to index 19
yhat=pamr$yhat[3] #choose threshold 6.763
str(pamr$prob) #it's a 3d dim array [i,j,k]=[nobvs,class,thresholdindex]
pprob=pamr$prob[,,3]
ytrue=SRBCT$y

#producing the output
vcrpamr=vcr.pamr.train(data=SRBCT, pamrfit=pamr, threshold= 6.7) #data is feeded in same format that pam accepts
vcrpamrcv=vcr.pamr.train(data=SRBCT, pamrfit=pamr, pamrfitcv=cvpamr, threshold= 8)
pamr.confusion(pamr, threshold=8)
vcrpamr$PAC


#silhouette visual plot
library(classmap)
?silplot #takes in a vcr out
silplot(vcrpamrcv, classLabels = c("EWS","BL","NB","RMS") ) #classLabels = c("EWS","BL","NB","RMS") for SRBCT
pamr.confusion(cvpamr, threshold=8)  # check if there is match with class label (about line of ordering label)
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
