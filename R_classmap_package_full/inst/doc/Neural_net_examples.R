## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
 fig.width  = 5 ,
 fig.height = 3.5,
 fig.align  = 'center'
)
oldpar <- list(mar = par()$mar, mfrow = par()$mfrow)

## -----------------------------------------------------------------------------
library(nnet)
library(classmap)

## -----------------------------------------------------------------------------
set.seed(123)
nn.iris <- nnet(Species ~ ., data = iris, size = 3)
names(nn.iris)
head(nn.iris$fitted.values) # matrix of posterior probabilities.
nn.iris$wts # the "weights" (coefficients, some are negative).
summary(nn.iris) 

## -----------------------------------------------------------------------------
train_x <- as.matrix(iris[, -5]) # original input: 4 variables.
# extract "weights" (coefficients) from each layer:
bias1  <- nn.iris$wts[c(1, 6, 11)] # intercepts to hidden layer
betas1 <- cbind(nn.iris$wts[2:5], nn.iris$wts[7:10], 
                nn.iris$wts[12:15]) # slopes to hidden layer
bias2  <- nn.iris$wts[c(16, 20, 24)] # intercepts to output layer
betas2 <- cbind(nn.iris$wts[17:19], nn.iris$wts[21:23], 
                nn.iris$wts[25:27]) # slopes to output layer.

## -----------------------------------------------------------------------------
sigmoid <- function(x) 1 / (1 + exp(-x)) 

H <- t(t(train_x %*% betas1) + bias1) # = hidden layer
X <- t(t(sigmoid(H) %*% betas2) + bias2) # = layer before softmax.
pairs(X, col = unclass(iris$Species)) # G=3 classes ==> 3 dimensions.

## -----------------------------------------------------------------------------
softmax <- function(x) exp(x) / sum(exp(x)) 
outprobs <- t(apply(X, 1, softmax))
trainprobs <- nn.iris$fitted.values # this is after softmax
range(trainprobs - outprobs) # near 0, so we have a match.

## -----------------------------------------------------------------------------
y <- iris[, 5]
vcrtrain <- vcr.neural.train(X, y, trainprobs)

names(vcrtrain)
vcrtrain$predint[c(1:10, 51:60, 101:110)] # prediction as integer
vcrtrain$pred[c(1:10, 51:60, 101:110)]    # prediction as label
vcrtrain$altint[c(1:10, 51:60, 101:110)]  # alternative label as integer
vcrtrain$altlab[c(1:10, 51:60, 101:110)]  # alternative label

# Probability of Alternative Class (PAC) of each object:
vcrtrain$PAC[1:3] 
#
summary(vcrtrain$PAC)

# f(i, g) is the distance from case i to class g:
vcrtrain$fig[1:5, ] # for the first 5 objects:

# The farness of an object i is the f(i, g) to its own class: 
vcrtrain$farness[1:5]
#
summary(vcrtrain$farness)

# The "overall farness" of an object is defined as the 
# lowest f(i, g) it has to any class g (including its own):
summary(vcrtrain$ofarness)

confmat.vcr(vcrtrain)

## -----------------------------------------------------------------------------
cols <- c("red", "darkgreen", "blue")

# stacked mosaic plot:
stplot <- stackedplot(vcrtrain, classCols = cols,
                     main = "Stacked plot of nnet on iris data") 
stplot

# silhouette plot:
silplot(vcrtrain, classCols = cols, main =
          "Silhouette plot of nnet on iris data") 

# class maps:
classmap(vcrtrain, "setosa", classCols = cols)
# Very tight class (low PAC, no high farness).
classmap(vcrtrain, "versicolor", classCols = cols)
# Not so tight, one point is predicted as virginica.
classmap(vcrtrain, "virginica", classCols = cols)
# Also tight.

## -----------------------------------------------------------------------------
test_x <- train_x[1:100, ]
ynew <- y[1:100]
ynew[c(1:10, 51:60)] <- NA
pairs(train_x, col = as.numeric(y) + 1, pch = 19) # 3 colors
pairs(test_x, col = as.numeric(ynew) + 1, pch = 19) # only red and green

## -----------------------------------------------------------------------------
predprobs <- predict(nn.iris, test_x, type = "raw")
range(predprobs - trainprobs[1:100, ]) # perfect match

# Reconstruct this prediction:
Hnew <- t(t(test_x %*% betas1) + bias1) # hidden layer
Xnew <- t(t(sigmoid(Hnew) %*% betas2) + bias2) # layer before softmax
probsnew <- t(apply(Xnew, 1, softmax))
range(probsnew - predprobs) # ~0, so we have a match

## -----------------------------------------------------------------------------
vcrtest <- vcr.neural.newdata(Xnew, ynew, probsnew, vcrtrain)

plot(vcrtest$predint, vcrtrain$predint[1:100]); abline(0, 1) # identical, OK
plot(vcrtest$altint, vcrtrain$altint[1:100]); abline(0, 1) # identical when not NA, OK
plot(vcrtest$PAC, vcrtrain$PAC[1:100]); abline(0, 1) # OK
vcrtest$farness # length 100, with NA's where ynew is NA
plot(vcrtest$farness, vcrtrain$farness[1:100]); abline(0, 1) # identical where not NA
plot(vcrtest$fig, vcrtrain$fig[1:100, ]); abline(0, 1) # same, OK
vcrtest$ofarness # 100, withOUT NA's: ofarness even exists for cases with missing ynew, for which farness cannot exist.
plot(vcrtest$ofarness, vcrtrain$ofarness[1:100]); abline(0, 1) # same, OK

confmat.vcr(vcrtest) # as expected:

## -----------------------------------------------------------------------------
stplot # to compare with:
stackedplot(vcrtest, classCols = cols, separSize = 1.5,
            minSize = 1, main = "Stacked plot of nnet on subset of iris data")
silplot(vcrtest, classCols = cols, main = "Silhouette plot of nnet on subset of iris data")

classmap(vcrtrain, 1, classCols = cols)
classmap(vcrtest, 1, classCols = cols) # same, but fewer points
#
classmap(vcrtrain, 2, classCols = cols)
classmap(vcrtest, 2, classCols = cols) # same, but fewer points
#
# classmap(vcrtrain, 3, classCols = cols)
# classmap(vcrtest, 3, classCols = cols)
# # Class number 3 with label virginica has no objects to visualize.

## -----------------------------------------------------------------------------
data("data_floralbuds")

## -----------------------------------------------------------------------------
set.seed(123)
nn.buds <- nnet(y ~ ., data = data_floralbuds, size = 2)
names(nn.buds)
head(nn.buds$fitted.values) # matrix of posterior probabilities of each class
nn.buds$wts # the "weights" (coefficients, some negative)
summary(nn.buds) # A hidden layer with 2 neurons:

## -----------------------------------------------------------------------------
train_x <- as.matrix(data_floralbuds[, -7])
# extract weights for each layer:
bias1  <- nn.buds$wts[c(1, 8)]
betas1 <- cbind(nn.buds$wts[2:7], nn.buds$wts[9:14])
bias2  <- nn.buds$wts[c(15, 18, 21, 24)]
betas2 <- cbind(nn.buds$wts[16:17], nn.buds$wts[19:20], 
                nn.buds$wts[22:23], nn.buds$wts[25:26])

## -----------------------------------------------------------------------------
H <- t(t(train_x %*% betas1) + bias1) # hidden layer
X <- t(t(sigmoid(H) %*% betas2) + bias2) # layer before softmax
outprobs <- t(apply(X, 1, softmax))
trainprobs <- nn.buds$fitted.values # posterior probabilities
range(trainprobs - outprobs) # ~0
# OK, we have reconstructed the posterior probabilities
pairs(X, col = data_floralbuds$y)

## -----------------------------------------------------------------------------
y <- data_floralbuds$y
vcrtrain <- vcr.neural.train(X, y, trainprobs)
cols <- c("saddlebrown", "orange", "olivedrab4", "royalblue3")

## -----------------------------------------------------------------------------
stackedplot(vcrtrain, classCols = cols, main =
              "Stacked plot of nnet on floral buds data")

# Silhouette plot:
silplot(vcrtrain, classCols = cols, 
        main = "Silhouette plot of nnet on floral buds data")



# Quasi residual plot:

PAC <- vcrtrain$PAC
feat <- rowSums(train_x); xlab <- "rowSums(X)"
# pdf("Floralbuds_quasi_residual_plot.pdf", width = 5, height = 4.8)
qresplot(PAC, feat, xlab = xlab, plotErrorBars = TRUE, fac = 2, main = "Floral buds: quasi residual plot")
# images with higher sum are easier to classify
# dev.off()
cor.test(feat, PAC, method = "spearman") 
# rho = -0.238, p-value < 2e-8 
# ==> decreasing trend is significant



# Class maps:

classmap(vcrtrain, "branch", classCols = cols)

classmap(vcrtrain, "bud", classCols = cols)

classmap(vcrtrain, "scales", classCols = cols)

classmap(vcrtrain, "support", classCols = cols)

