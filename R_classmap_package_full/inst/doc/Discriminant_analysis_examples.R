## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
 fig.width  = 5 ,
 fig.height = 3.5,
 fig.align  = 'center'
)
oldpar <- list(mar = par()$mar, mfrow = par()$mfrow)

## -----------------------------------------------------------------------------
library("classmap")
library("ggplot2")
library("gridExtra")

## -----------------------------------------------------------------------------
data(iris)
X <- iris[, 1:4]
y <- iris[, 5]
is.factor(y)
table(y)

pairs(X, col = as.numeric(y) + 1, pch = 19)

## -----------------------------------------------------------------------------
vcr.train <- vcr.da.train(X, y, rule = "QDA")
names(vcr.train)

## -----------------------------------------------------------------------------
vcr.train$predint 
vcr.train$pred[c(1:10, 51:60, 101:110)]
vcr.train$altint  
vcr.train$altlab[c(1:10, 51:60, 101:110)]

## -----------------------------------------------------------------------------
vcr.train$PAC[1:3] 
summary(vcr.train$PAC)

## -----------------------------------------------------------------------------
vcr.train$fig[1:5, ]

## -----------------------------------------------------------------------------
vcr.train$farness[1:5]
summary(vcr.train$farness)

## -----------------------------------------------------------------------------
summary(vcr.train$ofarness)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train, cutoff = 0.98)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train, showClassNumbers = TRUE)

## -----------------------------------------------------------------------------
cols <- c("red", "darkgreen", "blue")
stackedplot(vcr.train, classCols = cols, separSize = 1.5,
            minSize = 1, showLegend = TRUE)

stackedplot(vcr.train, classCols = cols, separSize = 1.5,
            minSize = 1, showLegend = TRUE, cutoff = 0.98)

## -----------------------------------------------------------------------------
stplot <- stackedplot(vcr.train, classCols = cols, 
                     separSize = 1.5, minSize = 1,
                     main = "QDA on iris data")
stplot

## -----------------------------------------------------------------------------
# pdf("Iris_QDA_silhouettes.pdf", width=5.0, height=4.6)
silplot(vcr.train, classCols = cols, 
        main = "Silhouette plot of QDA on iris data")      
# dev.off()

## -----------------------------------------------------------------------------
classmap(vcr.train, 1, classCols = cols)
classmap(vcr.train, 2, classCols = cols)
classmap(vcr.train, 3, classCols = cols) # With the default cutoff no farness values stand out:

## -----------------------------------------------------------------------------
# With a lower cutoff:
classmap(vcr.train, 3, classCols = cols, cutoff = 0.98)
# Now one point is to the right of the vertical line.
# It also has a black border, meaning that it is flagged
# as an outlier, in the sense that its farness to _all_
# classes is above 0.98.

## -----------------------------------------------------------------------------
Xnew <- X[c(1:50, 101:150), ]
ynew <- y[c(1:50, 101:150)]
ynew[c(1:10, 51:60)] <- NA
pairs(X, col = as.numeric(y) + 1, pch = 19) # 3 colors
pairs(Xnew, col = as.numeric(ynew) + 1, pch = 19) # only red and blue

## -----------------------------------------------------------------------------
vcr.test <- vcr.da.newdata(Xnew, ynew, vcr.train)

## -----------------------------------------------------------------------------
plot(vcr.test$predint, vcr.train$predint[c(1:50, 101:150)]); abline(0, 1)
plot(vcr.test$altint, vcr.train$altint[c(1:50, 101:150)]); abline(0, 1)
plot(vcr.test$PAC, vcr.train$PAC[c(1:50, 101:150)]); abline(0, 1)
vcr.test$farness 
plot(vcr.test$farness, vcr.train$farness[c(1:50, 101:150)]); abline(0, 1)
plot(vcr.test$fig, vcr.train$fig[c(1:50, 101:150), ]); abline(0, 1)
vcr.test$ofarness 
plot(vcr.test$ofarness, vcr.train$ofarness[c(1:50, 101:150)]); abline(0, 1)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.test)
confmat.vcr(vcr.test, cutoff = 0.98)

## -----------------------------------------------------------------------------
stplot # to compare with:
stackedplot(vcr.test, classCols = cols, separSize = 1.5, minSize = 1)

## -----------------------------------------------------------------------------
#pdf("Iris_test_QDA_silhouettes.pdf", width=5.0, height=4.3)
silplot(vcr.test, classCols = cols, 
        main = "Silhouette plot of QDA on iris subset") 
#dev.off()

## -----------------------------------------------------------------------------
classmap(vcr.train, 1, classCols = cols)
classmap(vcr.test, 1, classCols = cols) 

## ---- error = TRUE------------------------------------------------------------
classmap(vcr.train, 2, classCols = cols)
classmap(vcr.test, 2, classCols = cols)

## -----------------------------------------------------------------------------
classmap(vcr.train, 3, classCols = cols)
classmap(vcr.test, 3, classCols = cols) 

## -----------------------------------------------------------------------------
data(data_floralbuds)
X <- as.matrix(data_floralbuds[, 1:6])
y <- data_floralbuds$y
dim(X) # 550  6
length(y) # 550
table(y)
# branch     bud  scales support 
#     49     363      94      44 

# Pairs plot
cols <- c("saddlebrown", "orange", "olivedrab4", "royalblue3")
pairs(X, gap = 0, col = cols[as.numeric(y)]) # hard to separate visually

## -----------------------------------------------------------------------------
vcr.obj <- vcr.da.train(X, y)

## -----------------------------------------------------------------------------
confmat <- confmat.vcr(vcr.obj, showOutliers = FALSE)

confmat.vcr(vcr.obj) 

## -----------------------------------------------------------------------------
stackedplot(vcr.obj, classCols = cols, separSize = 0.6,
            minSize = 1.5,  main = "stacked plot of QDA on floral buds")

# Version in paper:
# pdf("Floralbuds_QDA_stackplot_without_outliers.pdf",
#     width=5, height=4.3)
# stackedplot(vcr.obj, classCols = cols, separSize = 0.6,
#             minSize = 1.5, showOutliers = FALSE,
#             htitle = "given class", vtitle = "predicted class")
# dev.off()

## -----------------------------------------------------------------------------
#pdf("Floralbuds_QDA_silhouettes.pdf", width=5.0, height=4.3)
silplot(vcr.obj, classCols = cols,
        main = "Silhouette plot of QDA on floral bud data")      
#dev.off()

## -----------------------------------------------------------------------------
PAC <- vcr.obj$PAC
feat <- rowSums(X); xlab = "rowSums(X)"
# pdf("Floralbuds_QDA_quasi_residual_plot.pdf", width=5, height=4.8)
qresplot(PAC, feat, xlab = xlab, plotErrorBars = TRUE, fac = 2, 
         main = "Floral buds: quasi residual plot")
# dev.off()

cor.test(feat, PAC, method = "spearman") 

## -----------------------------------------------------------------------------
labels <- c("branch", "bud", "scale", "support")

# classmap of class "bud"
#
# To identify the points that stand out:
# classmap(vcr.obj, 2, classCols = cols, identify = TRUE)
# Press "Esc" to get out.
#
# pdf("Floralbuds_QDA_classmap_bud.pdf", width=7, height=7)
par(mar = c(3.6, 3.5, 2.4, 3.5))
coords <- classmap(vcr.obj, 2, classCols = cols,
         main = "predictions of buds",
         cex = 1.5, cex.lab = 1.5, cex.axis = 1.5,
         cex.main = 1.5) 
# For marking points:
indstomark <- c(294, 70, 69, 152, 204) # from identify = TRUE above
labs  <- letters[seq_len(5)]
xvals <- coords[indstomark, 1] +
  c(0, 0.10, 0.14, 0.10, 0.08) # visual finetuning
yvals <- coords[indstomark, 2] +
  c(0.04, 0.04, 0, -0.03, +0.04)
text(x = xvals, y = yvals, labels = labs, cex = 1.5)
legend("topleft", fill = cols[1:4], legend = labels, 
       cex = 1, ncol = 1, bg = "white")
# dev.off()
par(oldpar)

## ---- fig.height = 8, fig.width = 8-------------------------------------------

#
# pdf(file = "Floralbuds_all_class_maps.pdf", width = 7, height = 7)
par(mfrow = c(2, 2))
par(mar = c(3.3, 3.2, 2.7, 1.0))
classmap(vcr.obj, 1, classCols = cols,
         main = "predictions of branches")
legend("topright", fill = cols, legend = labels,
       cex = 1, ncol = 1, bg = "white")
#
par(mar = c(3.3, 0.5, 2.7, 0.3))
classmap(vcr.obj, 2, classCols = cols,
         main = "predictions of buds")
labs  <- letters[seq_len(5)]
xvals <- coords[indstomark, 1] +
  c(0, 0.10, 0.14, 0.10, 0.08) # visual finetuning
yvals <- coords[indstomark, 2] +
  c(0.04, 0.04, 0, -0.03, 0.04)
# xvals <- c( 1.75, 1.68, 1.25, 3.25, 4.00)
# yvals <- c(0.045, 0.92, 0.54, 0.97, 0.045)
text(x = xvals, y = yvals, labels = labs, cex = 1.0)
legend("topleft", fill = cols, legend = labels,
       cex = 1, ncol = 1, bg = "white")
#
par(mar = c(3.3, 3.2, 2.7, 1.0))
classmap(vcr.obj, 3, classCols = cols,
         main = "predictions of scales")
legend("left", fill = cols, legend = labels,
       cex = 1, ncol = 1, bg = "white")
# 
par(mar = c(3.3, 0.5, 2.7, 0.3))
classmap(vcr.obj, 4, classCols = cols,
         main = "predictions of supports")
legend("topright", fill = cols, legend = labels,
       cex = 1, ncol = 1, bg = "white")
# dev.off()
par(oldpar)

## -----------------------------------------------------------------------------

mnist_url <- "https://wis.kuleuven.be/statdatascience/robust/data/mnist-rdata"
url.exists <- suppressWarnings(try(open.connection(url(mnist_url), open = "rt", timeout = 2),  silent = TRUE)[1], classes = "warning")

if (is.null(url.exists)) {load(url(mnist_url))} else {
  print(paste("The data source ", mnist_url, "is not active at the moment. The example can nevertheless be reproduced by downloading the mnist data from another source, formatting the training data to dimensions 60000 x 28 x 28, and running the code below."))
}
close(url(mnist_url))
X_train <- mnist$train$x
y_train <- as.factor(mnist$train$y)

head(y_train)
# Levels: 0 1 2 3 4 5 6 7 8 9
dim(X_train) # 60000    28    28
length(y_train) # 60000

## ---- fig.height = 1, fig.width = 1-------------------------------------------

plotImage = function(tempImage) {
  tdm = reshape2::melt(apply((tempImage), 2, rev))
  p = ggplot(tdm, aes(x = Var2, y = Var1, fill = (value))) +
    geom_raster() +
    guides(color = "none", size = "none", fill = "none") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_fill_gradient(low = "white", high = "black")
  p
}

plotImage(X_train[1, , ])
plotImage(X_train[2, , ])
plotImage(X_train[3, , ])

## ---- fig.height = 2, fig.width = 5-------------------------------------------
# Change the dimensions of X for the sequel:
dim(X_train) <- c(60000, 28 * 28)
dim(X_train) # 60000    784

# Sampled digit images:
set.seed(123)
sampledigits <- list()
for (i in 0:9) {
  digit <- i
  idx <- sample(which(y_train == digit), size = 1)
  tempImage <- matrix(unlist(X_train[idx, ]), 28, 28)
  sampledigits[[i + 1]] <- plotImage(tempImage) 
}
psampledigits <- grid.arrange(grobs = sampledigits, ncol = 5)
# ggsave("MNIST_sampled_images.pdf", plot = psampledigits,
#        width = 10, height = 1)


# Averaged digit images:
meanPlots <- list()
for (j in 0:9) {
  m.out <- colMeans(X_train[which(y_train == j), ])
  dim(m.out) <- c(28, 28)
  meanPlots[[j + 1]] <- plotImage(m.out) 
}
meanplot <- grid.arrange(grobs = meanPlots, ncol = 5)
# ggsave("MNIST_averaged_images.pdf", plot = meanplot,
#        width = 10, height = 1)

## -----------------------------------------------------------------------------
library(svd)
ptm <- proc.time()
svd.out <- svd::propack.svd(X_train, neig = 50)
(proc.time() - ptm)[3]
loadings <- svd.out$v
rm(svd.out)
dataProj <- as.matrix(X_train %*% loadings)
dim(dataProj)

## -----------------------------------------------------------------------------
vcr.train <- vcr.da.train(X = dataProj, y_train)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.train, showOutliers = FALSE)

cols <- c("red3", "darkorange", "gold2", "darkolivegreen3",
         "darkolivegreen4", "cadetblue3", "deepskyblue4", 
         "darkslateblue", "darkorchid3", "deeppink4")

# stacked plot in paper:
# pdf("MNIST_stackplot_with_outliers.pdf", width=5, height=4.3)
stackedplot(vcr.train, classCols = cols, separSize = 0.6,
            minSize = 1.5, htitle = "given class",
            main = "Stacked plot of QDA on MNIST training data", vtitle = "predicted class")
# dev.off()

## -----------------------------------------------------------------------------
# pdf("MNIST_QDA_silhouettes.pdf", width=5.0, height=4.6)
silplot(vcr.train, classCols = cols,
        main = "Silhouette plot of QDA on MNIST training data")      
# dev.off()

## -----------------------------------------------------------------------------
wnq <- function(string, qwrite=TRUE) { # auxiliary function
  # writes a line without quotes
  if (qwrite) write(noquote(string), file = "", ncolumns = 100)
}

showdigit <- function(digit=digit, i, plotIt = TRUE) {
  idx = which(y_train == digit)[i]
  # wnq(paste("Estimated digit: ", as.numeric(vcr.train$pred[idx]), sep=""))
  tempImage <- matrix(unlist(X_train[idx, ]), 28, 28)
  if (plotIt) {plot(plotImage(tempImage))}
  return(plotImage(tempImage))
}


## ---- fig.height = 4, fig.width = 6-------------------------------------------
digit <- 0
#
# To identify outliers:
# classmap(vcr.train, digit+1, classCols = cols, identify = TRUE)
# Press "Esc" to get out.
#
# pdf(paste0("MNIST_classmap_digit", digit, ".pdf"), width = 7, height = 7)
par(mar = c(3.6, 3.5, 2.4, 3.5))
coords <- classmap(vcr.train, digit + 1, classCols = cols,
         main = paste0("predictions of digit ",digit),
         cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
         cex.main = 1.5)
indstomark <- c(4000, 3964, 5891, 2485, 822, 
               2280, 2504, 3906, 5869, 1034) # from identify = TRUE
labs  <- letters[1:length(indstomark)]
xvals <- coords[indstomark, 1] +
  c(-0.04, -0.01, 0, -0.11, 0.06,
    0.07, 0.06, 0.10, 0.06, 0.09)
yvals <- coords[indstomark, 2] +
  c(-0.03, -0.03, -0.03, 0.022, -0.025, 
    -0.025, -0.035, -0.025, 0.03, 0.03)
text(x = xvals, y = yvals, labels = labs, cex = 1.5)
legend("left", fill = cols,
       legend = 0:9, cex = 1, ncol = 2, bg = "white")
# dev.off()
par(oldpar)

## ---- fig.height = 2, fig.width = 5-------------------------------------------
pred <- vcr.train$pred # needed for discussion plots
tempPreds <- (pred[which(y_train == digit)])[indstomark]
discussionPlots <- list()
for (i in 1:length(indstomark)) {
  idx <- indstomark[i]
  tempplot <- showdigit(digit, idx, plotIt = FALSE)
  tempplot <- arrangeGrob(tempplot, 
    bottom = paste0("(", labs[i], ") \"", tempPreds[i], "\""))
  discussionPlots[[i]] = tempplot
}
discussionPlot <- grid.arrange(grobs = discussionPlots, ncol = 5)
# ggsave(paste0("MNIST_discussionplot_digit", digit, ".pdf"),
#        plot = discussionPlot, width = 5,
#        height = (length(indstomark) %/% 5 +
#                    (length(indstomark) %% 5 > 0)))

## ---- fig.height = 4, fig.width = 6-------------------------------------------

digit <- 1
# pdf(paste0("MNIST_classmap_digit", digit, ".pdf"), width = 7, height = 7)
par(mar = c(3.6, 3.5, 2.4, 3.5))
classmap(vcr.train, digit + 1, classCols = cols,
         main = paste0("predictions of digit ", digit),
         cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
         cex.main = 1.5)
legend("left", fill = cols,
       legend = 0:9, cex = 1, ncol = 2, bg = "white")
# dev.off()
par(oldpar)

## ---- fig.height = 13, fig.width = 8------------------------------------------
# indices of the 1s predicted as 2 (takes a while):
#
indstomark <- which(vcr.train$predint[which(y_train == digit)] == 3)
length(indstomark) # 104
labs  <- letters[1:length(indstomark)]
pred <- vcr.train$pred # needed for discussion plots
tempPreds    <- (pred[which(y_train == digit)])[indstomark]
discussionPlots <- list()
for (i in 1:length(indstomark)) {
  idx <- indstomark[i]
  tempplot <- showdigit(digit, idx, FALSE)
  tempplot <- arrangeGrob(tempplot, 
        bottom = paste0("\"", tempPreds[i], "\""))
  discussionPlots[[i]] <- tempplot
}
discussionPlot <- grid.arrange(grobs = discussionPlots, 
                                         ncol = 8)
# ggsave(paste0("MNIST_discussionplot_digit", digit, "predictedAs2b.pdf"),
#        plot = discussionPlot, width = 10,
#        height = (length(indstomark) %/% 10 +
#                    (length(indstomark) %% 10 > 0)))
                 
# The digits 1 predicted as a 2 are mostly ones written with
# a horizontal line at the bottom.

## ---- fig.height = 4, fig.width = 6-------------------------------------------
digit <- 2
# To identify outliers:
# classmap(vcr.train, digit + 1, classCols = cols, identify = TRUE)
# Press "Esc" to get out.
#
# pdf(paste0("MNIST_classmap_digit", digit,".pdf"), width = 7, height = 7)
par(mar = c(3.6, 3.5, 2.4, 3.5))
coords <- classmap(vcr.train, digit + 1, classCols = cols,
                  main = paste0("predictions of digit", digit), cex = 1.5, cex.lab = 1.5, cex.axis = 1.5,
                  cex.main = 1.5)
indstomark <- c(3164, 5434, 2319 , 4224, 3682, 
               2642, 4920, 1233, 3741, 3993) # from identify = TRUE
labs  <- letters[1:length(indstomark)]
xvals <- coords[indstomark, 1] +
  c(0, 0.08, 0, 0, 0, 0, 0, 0, 0, 0)
yvals <- coords[indstomark, 2] +
  c(-0.03, -0.03, -0.03, -0.03, -0.03, 
    -0.03, -0.03, -0.03, 0.03, 0.03)  
text(x = xvals, y = yvals, labels = labs, cex = 1.5)
legend("right", fill = cols,
       legend = 0:9, cex = 1, ncol = 2, bg = "white")
# dev.off()
par(oldpar)

## ---- fig.height = 2, fig.width = 5-------------------------------------------
pred <- vcr.train$pred # needed for discussion plots
tempPreds    <- (pred[which(y_train == digit)])[indstomark]
discussionPlots <- list()
for (i in 1:length(indstomark)) {
  idx <- indstomark[i]
  tempplot <- showdigit(digit, idx, FALSE)
  tempplot <- arrangeGrob(tempplot, 
        bottom = paste0("(", labs[i], ") \"", tempPreds[i], "\""))
  discussionPlots[[i]] <- tempplot
}
discussionPlot <- grid.arrange(grobs = discussionPlots, 
                                         ncol = 5)
# ggsave(paste0("MNIST_discussionplot_digit", digit, ".pdf"),
#        plot = discussionPlot, width = 5,
#        height = (length(indstomark) %/% 5 +
#                    (length(indstomark) %% 5 > 0)))

## ---- fig.height = 1, fig.width = 1-------------------------------------------
X_test <- mnist$test$x
y_test <- as.factor(mnist$test$y)

head(y_test)
#
dim(X_test) # 10000    28    28
length(y_test) # 10000

plotImage(X_test[1, , ])
plotImage(X_test[2, , ])
plotImage(X_test[3, , ])

dim(X_test) <- c(10000, 28 * 28)
dim(X_test) # 10000  784

dataProj_test <- as.matrix(X_test %*% loadings)

## -----------------------------------------------------------------------------
vcr.test <- vcr.da.newdata(Xnew = dataProj_test,
                           ynew = y_test,
                           vcr.da.train.out = vcr.train)

## -----------------------------------------------------------------------------
confmat.vcr(vcr.test, showOutliers = FALSE, showClassNumbers = TRUE)

# In supplementary material:
# pdf("MNISTtest_stackplot_with_outliers.pdf", width = 5, height = 4.3)
stackedplot(vcr.test, classCols = cols, separSize = 0.6,
            main = "Stacked plot of QDA on MNIST test data",
            minSize = 1.5)
# dev.off()

## -----------------------------------------------------------------------------
#pdf("MNIST_test_QDA_silhouettes.pdf", width = 5.0, height = 4.6)
silplot(vcr.test, classCols = cols,
        main = "Silhouette plot of QDA on MNIST test data")      
#dev.off()

## ---- fig.height = 4, fig.width = 6-------------------------------------------
showdigit_test <- function(digit = digit, i, plotIt = TRUE) {
  idx = which(y_test == digit)[i]
  # wnq(paste("Estimated digit: ", as.numeric(vcr.test$pred[idx]), sep = ""))
  tempImage <- matrix(unlist(X_test[idx, ]), 28, 28)
  if (plotIt) {plot(plotImage(tempImage))}
  return(plotImage(tempImage))
}

digit <- 0
# classmap(vcr.test, digit+1, classCols = cols, identify = TRUE)
# pdf(paste0("MNISTtest_classmap_digit", digit,".pdf"))
par(mar = c(3.6, 3.5, 2.4, 3.5))
coords <- classmap(vcr.test, digit + 1, classCols = cols,
         main = paste0("predictions of digit ", digit),
         cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
         cex.main = 1.5)
indstomark <- c(140, 630, 241, 967, 189,
               377, 78, 943, 64, 354)
labs  <- letters[1:length(indstomark)]
xvals <- coords[indstomark, 1] +
  c(0.08, 0.07, -0.07, 0.06, 0,
    0.04, 0.05, 0.09, -0.04, 0.09)
yvals <- coords[indstomark, 2] +
  c(-0.025, -0.03, -0.024, -0.025, -0.03, 
    -0.03, -0.03, 0.022, 0.035, 0.03)
text(x = xvals, y = yvals, labels = labs, cex = 1.5)
legend("left", fill = cols,
       legend = 0:9, cex = 1, ncol = 2, bg = "white")
# dev.off()
par(oldpar)

## ---- fig.height = 2, fig.width = 5-------------------------------------------
pred <- vcr.test$pred # needed for discussion plots
tempPreds <- (pred[which(y_test == digit)])[indstomark]
discussionPlots <- list()
for (i in 1:length(indstomark)) {
  idx <- indstomark[i]
  tempplot <- showdigit_test(digit, idx, FALSE)
  tempplot <- arrangeGrob(tempplot, 
      bottom = paste0("(", labs[i], ") \"", tempPreds[i], "\""))
  discussionPlots[[i]] <- tempplot
}
discussionPlot <- grid.arrange(grobs = discussionPlots, 
                                         ncol = 5)
# ggsave(paste0("MNISTtest_discussionplot_digit", digit, ".pdf"),
#        plot = discussionPlot, width = 5,
#        height = (length(indstomark) %/% 5 +
#                    (length(indstomark) %% 5 > 0)))

## ---- fig.height = 4, fig.width = 6-------------------------------------------
digit <- 3
# classmap(vcr.test, digit + 1, classCols = cols, identify = TRUE)
# pdf(paste0("MNISTtest_classmap_digit", digit, ".pdf"))
par(mar = c(3.6, 3.5, 2.4, 3.5))
coords <- classmap(vcr.test, digit + 1, classCols = cols,
         main = paste0("predictions of digit ", digit),
         cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
         cex.main = 1.5)
indstomark <- c(883, 659, 262, 60, 310,
               832, 223, 784, 835, 289)
labs  <- letters[1:length(indstomark)]
xvals <- coords[indstomark, 1] +
  c(-0.01, 0.08, -0.10, 0.06, 0.07, 
    0.06, 0.03, 0.11, 0.02, 0.06)
yvals <- coords[indstomark, 2] +
  c(0.035, 0.033, -0.017, -0.022, -0.025, 
    -0.025, -0.033, -0.022, 0.035, 0.038)
text(x = xvals, y = yvals, labels = labs, cex = 1.5)
legend("right", fill = cols,
       legend = 0:9, cex = 1, ncol = 2, bg = "white")
# dev.off()
par(oldpar)

## ---- fig.height = 2, fig.width = 5-------------------------------------------
pred <- vcr.test$pred # needed for discussion plots
tempPreds    <- (pred[which(y_test == digit)])[indstomark]
discussionPlots <- list()
for (i in 1:length(indstomark)) {
  idx <- indstomark[i]
  tempplot <- showdigit_test(digit, idx, FALSE)
  tempplot <- arrangeGrob(tempplot, 
    bottom = paste0("(", labs[i], ") \"", tempPreds[i], "\""))
  discussionPlots[[i]] <- tempplot
}
discussionPlot <- grid.arrange(grobs = discussionPlots, 
                                         ncol = 5)
# ggsave(paste0("MNISTtest_discussionplot_digit", digit, ".pdf"),
#        plot = discussionPlot, width = 5,
#        height = (length(indstomark) %/% 5 +
#                    (length(indstomark) %% 5 > 0)))

