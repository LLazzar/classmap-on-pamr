## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
 fig.width  = 5 ,
 fig.height = 3.5,
 fig.align  = 'center'
)
oldpar <- list(mar = par()$mar, mfrow = par()$mfrow)

## -----------------------------------------------------------------------------
library(randomForest)
library(classmap)

## -----------------------------------------------------------------------------
data("data_instagram")
traindata <- data_instagram[which(data_instagram$dataType == "train"), -13]
str(traindata)


# The variable names and their interpretation are
colnames(traindata)

# profile.pic: binary, indicates whether profile has picture
# nums.length.username: ratio of number of numerical chars in username to its length
# fullname.words: number of words in full name
# nums.length.fullname: ratio of number of numerical characters in full name to its length
# name..username: binary, indicates whether name == username of the profile
# description.length: length of the description/biography of the profile (in number of characters)
# external.URL: binary, indicates whether profile has external url
# private: binary, indicates whether profile is private or not
# X.posts: number of posts made by profile
# X.followers: number of followers
# X.follows: numbers of follows
# y: whether profile is fake or not.

x_train <- traindata[, -12]
y_train <- traindata[, 12]

dim(traindata)
table(traindata$y) # 50/50 split of genuine/fake accounts:

## -----------------------------------------------------------------------------
set.seed(71) 
rfout <- randomForest(y ~ ., data = traindata, keep.forest = TRUE)

## -----------------------------------------------------------------------------
mytype <- list(symm = c(1, 5, 7, 8)) 

## -----------------------------------------------------------------------------
vcrtrain <- vcr.forest.train(X = x_train, y = y_train,
                            trainfit = rfout, type = mytype)

names(vcrtrain)
vcrtrain$predint[c(1:10, 301:310)] # prediction as integer
vcrtrain$pred[c(1:10, 301:310)]    # prediction as label
vcrtrain$altint[c(1:10, 301:310)]  # alternative label as integer
vcrtrain$altlab[c(1:10, 301:310)]  # alternative label

# Probability of Alternative Class (PAC) of each object:
vcrtrain$PAC[1:3] 
#
summary(vcrtrain$PAC)

# f(i, g) is the distance from case i to class g:
vcrtrain$fig[1:3, ] # for the first 3 objects:

# The farness of an object i is the f(i, g) to its own class: 
vcrtrain$farness[1:3]
#
summary(vcrtrain$farness)

# The "overall farness" of an object is defined as the 
# lowest f(i, g) it has to any class g (including its own):
summary(vcrtrain$ofarness)

sum(vcrtrain$ofarness > 0.99, na.rm = TRUE) 
# With the default cutoff = 0.99 we find 6 outliers,
# also shown in the last column of the confusion matrix:

confmat.vcr(vcrtrain) 

# If we do not want to show the outliers:
confmat.vcr(vcrtrain, showOutliers = FALSE)

# Note that the accuracy is computed before any objects
# are flagged, so it does not depend on the cutoff.
# Here the accuracy is `perfect' due to overfitting. 
# The out-of-box prediction accuracy is about 92%.
cols <- c("blue", "red3")

## -----------------------------------------------------------------------------
stackedplot(vcrtrain, classCols = cols, main =
              "Instagram training data")

# Silhouette plot:
silplot(vcrtrain, classCols = cols)
# Here all the s(i) are nonnegative (due to overfitting).

# Class maps:
classmap(vcrtrain, "genuine", classCols = cols) #, identify = TRUE)

# farness outliers from furthest to closer: 45, 25, 41
x_train[c(45, 25, 41), ] # they have huge numbers of followers.

classmap(vcrtrain, "fake", classCols = cols) #, identify = TRUE)
# only case 261 is borderline far.

## -----------------------------------------------------------------------------
testdata <- data_instagram[which(data_instagram$dataType == "test"), -13]
Xnew <- testdata[, -12]
ynew <- testdata[, 12]

## -----------------------------------------------------------------------------
vcrtest <- vcr.forest.newdata(Xnew, ynew, vcrtrain)

confmat.vcr(vcrtest)

## -----------------------------------------------------------------------------
stackedplot(vcrtest, classCols = cols, 
            main = "RF on Instagram test data")

# Silhouette plot:
silplot(vcrtest, classCols = cols, main =
          "Silhouettes of RF on Instagram test data") # now some s(i) are negative

## -----------------------------------------------------------------------------

## Class of genuine accounts:

classmap(vcrtest, "genuine", classCols = cols) #, identify = TRUE)

# one farness outlier:
Xnew[c(30), ]
# has very lengthy bio/description
# has large number of X.posts
# has very large number of followers and follows

# genuine misclassified as fake: from highest PAC to lowest
Xnew[c(21, 29, 51), ] # and 2 more borderline cases
# They have some unusual characteristics for their class:
# * 21, 29 have a (very) high nums.length.username, i.e. the
#   percentage of numerical characters in the username.
# * 21, 29 have a full name of only 1 word.
# * 21, 29 and 51 have description.length = 0, i.e. no 
#   description/biography of their profile.
# * they all have low X.posts (even 0 for case 21), i.e.
#   relatively few previous posts.
# All of these characteristics are more common for fake profiles
# than for genuine profiles, as we can see below:

trcols <- cols[as.numeric(y_train)]

plot(x_train[, 1], col = trcols, main = "profile.pic")
# fakes are less likely to have a profile picture
plot(x_train[, 2], col = trcols, main = "nums.length.username")
# is higher for fakes
plot(x_train[, 3], col = trcols, main = "fullname.words")
# is lower for fakes
plot(x_train[, 4], col = trcols, main = "nums.length.fullname")
# is a bit higher for fakes
plot(x_train[, 5], col = trcols, main = "name..username")
# mostly 0 for genuine; fakes have a few values 1
plot(x_train[, 6], col = trcols, main = "description.length")
# fakes are typically lower, and more often zero
plot(x_train[, 7], col = trcols, main = "external.URL")
# fakes never had them, genuines sometimes did
plot(x_train[, 8], col = trcols, main = "private")
# no visible difference
plot((x_train[, 9])^0.1, col = trcols, main = "X.posts")
# fakes have fewer posts, and often none
plot((x_train[, 10])^0.1, col = trcols, main = "X.followers")
# fakes have fewer followers, sometimes none
plot((x_train[, 11])^0.1, col = trcols, main = "X.follows")
# fakes follow a bit fewer, but the difference is small.


## Class of fake accounts:

classmap(vcrtest, "fake", classCols = cols) #, identify = TRUE)

# Fake identified as genuine, from highest PAC to lower:
# c(27, 51, 34, 23, 58)
Xnew[which(ynew == "fake")[c(27, 34, 51, 23, 58)], ]

# These have a number of characteristics which are more common 
# for genuine profiles:
#
# all have profile pictures
# none have numerical characters in username
# none have numerical characters in fullname
# 27 has a lengthy bio description
# all have a relatively high number of followers
# all have a relatively high number of follows.

