## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
fig.width  = 5 ,
fig.height = 3.5,
fig.align  = 'center'
)
oldpar <- list(mar = par()$mar, mfrow = par()$mfrow)

## -----------------------------------------------------------------------------
library("classmap")

## -----------------------------------------------------------------------------
data(iris)
X <- iris[, 1:4]
y <- iris[, 5]
is.factor(y) 
table(y) 

X <- scale(X) # if you want to scale, best to do it yourself.
pairs(X, col = as.numeric(y) + 1, pch = 19)
dis <- dist(X, method = "euclidean")

## -----------------------------------------------------------------------------
vcr.train <- vcr.knn.train(dis, y, k = 5)

vcr.train <- vcr.knn.train(X, y, k = 5) # Gives the same result.

## -----------------------------------------------------------------------------
names(vcr.train)
vcr.train$predint 
vcr.train$pred[c(1:10, 51:60, 101:110)]
vcr.train$altint  

## -----------------------------------------------------------------------------
vcr.train$altlab  # has NA's for the same cases:

## -----------------------------------------------------------------------------
vcr.train$k 
vcr.train$ktrues[1:10]

## -----------------------------------------------------------------------------
vcr.train$PAC # length 150, all between 0 and 1, no NA's

## -----------------------------------------------------------------------------
vcr.train$fig[1:5, ]

## -----------------------------------------------------------------------------
vcr.train$farness[1:5]
summary(vcr.train$farness)

## -----------------------------------------------------------------------------
summary(vcr.train$ofarness)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train, showOutliers = FALSE)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train, cutoff = 0.98)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train, showClassNumbers = TRUE)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 1))
cols <- c("red", "darkgreen", "blue")
stackedplot(vcr.train, classCols = cols, separSize = 1.5,
            minSize = 1, showOutliers = FALSE, showLegend = TRUE,
            main = "Stacked plot of kNN on iris data")

stackedplot(vcr.train, classCols = cols, separSize = 1.5,
            minSize = 1, showLegend = TRUE,
            main = "Stacked plot of kNN on iris data")
par(oldpar)

## -----------------------------------------------------------------------------
stplot <- stackedplot(vcr.train, classCols = cols, 
                      separSize = 1.5, minSize = 1,
                      main = "Stacked plot of kNN on iris data")
stplot 

## -----------------------------------------------------------------------------
# pdf("Iris_kNN_silhouettes.pdf", width=5.0, height=4.6)
silplot(vcr.train, classCols = cols, 
        main = "Silhouette plot of kNN on iris data")      
# dev.off()

## -----------------------------------------------------------------------------
classmap(vcr.train, 1, classCols = cols)

## -----------------------------------------------------------------------------
classmap(vcr.train, 2, classCols = cols)

## -----------------------------------------------------------------------------
classmap(vcr.train, 3, classCols = cols)

## -----------------------------------------------------------------------------
Xnew <- X[c(1:50, 101:150),]
ynew <- y[c(1:50, 101:150)]
ynew[c(1:10, 51:60)] <- NA
pairs(X, col = as.numeric(y) + 1, pch = 19) # 3 colors
pairs(Xnew, col = as.numeric(ynew) + 1, pch = 19) # only red and blue

## -----------------------------------------------------------------------------
vcr.test <- vcr.knn.newdata(Xnew, ynew, vcr.train, LOO = TRUE)

## -----------------------------------------------------------------------------
plot(vcr.test$predint, vcr.train$predint[c(1:50, 101:150)]); abline(0, 1)
plot(vcr.test$altint, vcr.train$altint[c(1:50, 101:150)]); abline(0, 1)
plot(vcr.test$PAC, vcr.train$PAC[c(1:50, 101:150)]); abline(0, 1) 
vcr.test$farness 
plot(vcr.test$farness, vcr.train$farness[c(1:50, 101:150)]); abline(0, 1)
plot(vcr.test$fig, vcr.train$fig[c(1:50, 101:150), ]); abline(0, 1)
vcr.test$ofarness  # This does exist for every case, even if is given label is NA:
plot(vcr.test$ofarness, vcr.train$ofarness[c(1:50, 101:150)]); abline(0, 1)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.test) 

stplot 
stackedplot(vcr.test, classCols = cols, separSize = 1.5, minSize = 1, main = "Stacked plot of kNN on iris subset") 

## -----------------------------------------------------------------------------
# pdf("Iris_test_kNN_silhouettes.pdf", width=5.0, height=4.6)
silplot(vcr.test, classCols = cols, 
        main = "Silhouette plot of kNN on iris subset")      
# dev.off()

## -----------------------------------------------------------------------------
classmap(vcr.train, 1, classCols = cols)
classmap(vcr.test, 1, classCols = cols) 

## ---- error = TRUE------------------------------------------------------------
classmap(vcr.train, 2, classCols = cols)
classmap(vcr.test, 2, classCols = cols)
# Class number 2 with label versicolor has no objects to visualize.

## -----------------------------------------------------------------------------
classmap(vcr.train, 3, classCols = cols)
classmap(vcr.test, 3, classCols = cols) # same, but fewer points

## -----------------------------------------------------------------------------
library(kernlab)
#?kernlab::spam 
data(spam)

## -----------------------------------------------------------------------------
y <- spam$type
table(y)
X <- spam[, -58] 
sum(duplicated(X))

## -----------------------------------------------------------------------------
vcr.obj <- vcr.knn.train(scale(X), y, k = 5) 
names(vcr.obj) 
vcr.obj$predint[1:100]; length(vcr.obj$predint) 
vcr.obj$altint[1:100]
vcr.obj$k 
summary(vcr.obj$ktrues) 

## -----------------------------------------------------------------------------
confmat.vcr(vcr.obj, showOutliers = FALSE)

confmat.vcr(vcr.obj)

## -----------------------------------------------------------------------------
cols <- c("deepskyblue2", "red") 
# nonspam is blue, spam is red

stackedplot(vcr.obj, separSize = 1.5, minSize = 1,
            classCols = cols, showOutliers = FALSE,
            main = "Spam data")

stackedplot(vcr.obj, separSize = 1.5, minSize = 2,
            classCols = cols, main = "Spam data")

## -----------------------------------------------------------------------------
# pdf("Spam_kNN_silhouettes.pdf", width=5.0, height=4.6)
silplot(vcr.obj, classCols = cols, 
        main = "Silhouette plot of kNN on spam data")      
# dev.off()

## ---- fig.height = 4, fig.width = 6-------------------------------------------
# To identify the points that stand out:
# classmap(vcr.obj, 1, classCols = cols, identify = T)
# Press "Esc" to get out.
#
# Class map in paper:
# pdf("Spamdata_classmap_ham.pdf")
par(mar = c(3.5, 4.3, 2.0, 0.2))
coords <- classmap(vcr.obj, 1, classCols = cols,
                   main = "predictions of non-spam mails",
                   cex = 1.5, cex.lab = 1.5, cex.axis = 1.5,
                   cex.main = 1.5, maxfactor = 1.03)
# From identify = TRUE above we can mark points:
indstomark <- c(202, 1434, 1596, 2651, 1576, 1804)
labs  <- letters[1:length(indstomark)]
xvals <- coords[indstomark, 1] +
  c(0, 0.125, 0.125, 0.125, 0.125,0.0625) # visual finetuning
yvals <- coords[indstomark, 2] +  
  c(-0.03, 0, 0, 0, 0, 0.03)
text(x = xvals, y = yvals, labels = labs, cex = 1.5)
legend("topleft", fill = cols,
       legend = c("ham", "spam"), cex = 1,
       ncol = 1, bg = "white")
# dev.off()
par(oldpar)
# To interpret the marked points:
# markInds = which(y == "nonspam")[indstomark]
# X[markInds, ] 

## ---- fig.height = 4, fig.width = 6-------------------------------------------
#
# To identify the points that stand out:
# classmap(vcr.obj, 2, classCols = cols, identify = TRUE)
# Press "Esc" to get out.
#
# Class map in paper:
# pdf("Spamdata_classmap_spam.pdf")
par(mar = c(3.5, 4.3, 2.0, 0.2))
coords <- classmap(vcr.obj, 2, classCols = cols,
                   main = "predictions of spam mails",
                   cex = 1.5, cex.lab = 1.5, cex.axis = 1.5,
                   cex.main = 1.5, maxfactor = 1.03)
indstomark <- c(1306, 1294, 1754, 177)
labs  <- letters[6 + (1:length(indstomark))]
xvals <- coords[indstomark, 1] + c(0.1, 0.1, 0.1, 0.1)
yvals <- coords[indstomark, 2] + c(-0.03, 0, 0, 0.03)
text(x = xvals, y = yvals, labels = labs, cex = 1.5)
legend("topleft", fill = cols,
       legend = c("ham", "spam"), cex = 1,
       ncol = 1, bg = "white")
# dev.off()
par(oldpar)

