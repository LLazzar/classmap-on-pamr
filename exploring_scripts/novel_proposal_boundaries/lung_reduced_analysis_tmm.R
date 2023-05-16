labels=read.csv("sample_datasets/GSE62944/lung_subset/labels.csv")
datar=read.csv("sample_datasets/GSE62944/lung_subset/FC_expression.csv")
covariates=read.csv("sample_datasets/GSE62944/lung_subset/covariates.csv")

## Perform TMM normalization ####
library(edgeR)
d <- DGEList(as.matrix(datar[,2:ncol(datar)]))
d <- calcNormFactors(d, method="TMM")
normCounts <- cpm(d)
###################

## Data Prep ########
data=list()
data$y=as.matrix(labels[,2])
data$x=normCounts
str(data$x)
str(data$y)
###############

library(Rtsne)
?Rtsne
tsne_results <- Rtsne(t(data$x), dims = 2, perplexity = 30)
plot(tsne_results$Y,col=col_vec_true, pch=21, bg=intens)
plot(tsne_results$Y,col=intens, pch=20, cex=2)

mds=cmdscale(vcrpamr$pwd, k=2)
mds=cmdscale(pairwiser, k=2)
plot(mds, col=col_vec_true, pch=21, bg=intens, cex=0.8) #how seen by algorithm by trues belonging
plot(mds, pch=20)


plot(mds, col=col_vec_true, pch=20, cex=0.8)
plot(mds, col=col_vec, pch=20, cex=0.8)

intens=rep(NA,length(col_vec_true))
for (i in 1:length(col_vec_true)){
  pal <- colorRamp(c(col_vec_true[i], "white"))
  intens[i] = rgb(pal(0.9*vcrpamr$PAC[i]), maxColorValue = 255)
}

pal <- colorRamp(c("blue", "white"))
pal(0.6)
rgb(pal(0.6),maxColorValue = 255)

#FITTING

library(pamr)

pamr=pamr.train(data)
thres=pamr.adaptthresh(pamr)
pamr=pamr.train(data, threshold.scale = thres)
pamr
pamrcv=pamr.cv(pamr, data)
pamrcv

labels=unlist(pamr$yhat[14], use.names = FALSE)
colors <- c("red", "green", "blue")
col_vec <- colors[match(labels, c("LUAD", "LUSC", "NORM"))]

labels=as.numeric(as.factor(labels))

col_vec_true = colors[match(data$y, c("LUAD", "LUSC", "NORM"))]




source("feature_code/R/VCR_pamr.R")

#VISUAL TOOLS

library(classmap)
vcrpamr=vcr.pamr.train(data,pamr, pamrfitcv=NULL, threshold=23)


#calculating pairwise distance all variables

sd=vcrpamr$sd
n=ncol(data$x)
pairwise=matrix(NA, nrow = n, ncol = n)
x<- t(apply(t(data$x), 1, function(x) x / sd))
pairwise=dist(x)


#alternative with mahalobis to check
sd=vcrpamr$sd
posid=rep(TRUE,length(sd))
pwd=pw_mdS2(data$x, sd, weight=posid)

#calculating pairwise distance variables selected
posid=vcrpamr$posid
sd=vcrpamr$sd
n=ncol(data$x)
pairwise=matrix(NA, nrow = n, ncol = n)
sdr=sd[posid]
xr=data$x[posid,]
pairwiser=matrix(NA, nrow = n, ncol = n)
xr<- t(apply(t(xr), 1, function(x) x / sdr))
pairwiser=dist(xr)
pairwiser= as.matrix(dist(xr)^2)




#######


#######

?silplot
silplot(vcrpamr)
classmap(vcrpamr, whichclass = 1)

#QUASI RESIDUAL

gender=unlist(covariates[40,2:ncol(covariates)], use.names = FALSE)
gender_coded <- ifelse(gender == "MALE", 0, 1)

#
average_expression_level=colSums(data$x)

#
age=as.numeric(unlist(covariates[2,2:ncol(covariates)], use.names = FALSE))

#
carbon_test=as.numeric(unlist(covariates[22,2:ncol(covariates)], use.names = FALSE))
#how are NAs treated in qresplot?? #test revealed that it ignores them all good!

#
organ_local=unlist(covariates[6,2:ncol(covariates)], use.names = FALSE)
lookup <- c("[Discrepancy]" = NA, "[Not Available]" = NA, "Bronchial" = 1, "L-Lower" = 2,
            "L-Upper" = 3, "Other (please specify)" = NA, "R-Lower" = 4, "R-Middle" = 5,
            "R-Upper" = 6)
coded_organ <- factor(organ_local, levels = names(lookup))
coded_organ <- as.numeric(lookup[coded_organ])

#
stage=unlist(covariates[33,2:ncol(covariates)], use.names = FALSE)
lookup <- c("[Discrepancy]" = NA, "[Not Available]" = NA, "Stage I" = 1, "Stage IA" = 2,
            "Stage IB" = 3, "Stage II" = 4, "Stage IIA" = 5, "Stage IIB" = 6, "Stage III" = 7,
            "Stage IIIA" = 8, "Stage IIIB" = 9, "Stage IV" = 10)
coded_stage <- factor(stage, levels = names(lookup))
coded_stage <- as.numeric(lookup[coded_stage])

#
detail=unlist(covariates[16,2:ncol(covariates)], use.names = FALSE)
table(detail)

#PLOTTING
?qresplot
dev.off()
qres=qresplot(vcrpamr$PAC, coded_stage, xlab = "Age (years)", opacity = 0.5,
              main = "quasi residual plot for male passengers", plotErrorBars = TRUE)


#testing missing values
nona=which(!is.na(carbon_test))
vcrpamr$PAC[nona]

dev.off()
qres=qresplot(vcrpamr$PAC[400:800], carbon_test[400:800], xlab = "Age (years)", opacity = 0.5,
              main = "quasi residual plot for male passengers", plotLoess = TRUE, identify = TRUE, cex=0.8)


#disturbance genes

coeff=rep(NA,ncol(t(data$x)))

for (i in 1:length(coeff)) {
  lm=lm(vcrpamr$PAC~t(data$x)[,i])
  coeff[i]=lm$coefficients[2]
}

sorted_indices <- order(coeff, decreasing = TRUE)
sorted_indices[1:10]
which.max(coeff)
lm(vcrpamr$PAC~t(khanc$x)[,424])

qres=qresplot(vcrpamr$PAC, t(data$x)[,], xlab = "Age (years)", opacity = 0.5,
              main = "quasi residual plot for male passengers", plotLoess = TRUE)
