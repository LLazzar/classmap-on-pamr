## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
 fig.width  = 5 ,
 fig.height = 3.5,
 fig.align  = 'center'
)
oldpar <- list(mar = par()$mar, mfrow = par()$mfrow)

## -----------------------------------------------------------------------------
library(rpart)
library(classmap)

## -----------------------------------------------------------------------------
data("data_titanic")

traindata <- data_titanic[which(data_titanic$dataType == "train"), -13]

dim(traindata) 
colnames(traindata)
# SibSp: number of siblings + spouses aboard
# Parch: number of parents + children aboard

str(traindata)
table(traindata$y)

## -----------------------------------------------------------------------------
set.seed(123) 
rpart.out <- rpart(y ~ Pclass + Sex + SibSp + 
                    Parch + Fare + Embarked,
                  data = traindata, method = 'class',
                  model = TRUE)

## plot the tree:
# pdf(file = "titanic_train_tree.pdf", width = 6, height = 6)
rpart.plot::rpart.plot(rpart.out, box.palette = "RdBu")
# dev.off()

## -----------------------------------------------------------------------------
rpart.out$variable.importance

## -----------------------------------------------------------------------------
mytype <- list(nominal = c("Name", "Sex", "Ticket", "Cabin", "Embarked"), ordratio = c("Pclass"))

## -----------------------------------------------------------------------------
x_train <- traindata[, -12]
y_train <- traindata[,  12]
vcrtrain <- vcr.rpart.train(x_train, y_train, rpart.out, mytype)
names(vcrtrain)

vcrtrain$predint[1:10] # prediction as integer
vcrtrain$pred[1:10]    # prediction as label
vcrtrain$altint[1:10]  # alternative label as integer
vcrtrain$altlab[1:10]  # alternative label

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
# No farness is considered outlying in these data.

confmat.vcr(vcrtrain) 

cols <- c("firebrick", "blue")

## -----------------------------------------------------------------------------
stackedplot(vcrtrain, classCols = cols,
            main = "Stacked plot of rpart on Titanic training data")

## -----------------------------------------------------------------------------
silplot(vcrtrain, classCols = cols, 
        main = "Silhouettes of rpart on Titanic training data")
# silplot.out <- silplot(vcrtrain, classCols = cols)
# ggsave("titanic_train_silhouettes.pdf", silplot.out,
#        width = 5, height = 5)

## -----------------------------------------------------------------------------
hist(x_train$Age)

# Quasi residual plot versus age, for males only:

# pdf("titanic_qrp_versus_age_males.pdf", width = 4.8, height = 4.8)
PAC <- vcrtrain$PAC[which(x_train$Sex == "male")]
feat <- x_train$Age[which(x_train$Sex == "male")]
qresplot(PAC, feat, xlab = "Age (years)", opacity = 0.5,
         main = "quasi residual plot for male passengers",
         plotLoess = TRUE)
text(x = 14, y = 0.60, "loess curve", col = "red", cex = 1)
# dev.off()

## -----------------------------------------------------------------------------
classmap(vcrtrain, "casualty", classCols = cols)
# classmap(vcrtrain, "casualty", classCols = cols, identify = TRUE)

# blue points top right: cases "a" and "b" in the paper
x_train[which(y_train == "casualty")[119], ]
# Woman in 1st class, should have survived.
#
x_train[which(y_train == "casualty")[192], ]

# Similar, but child.

# red point most to the right: case "c" in the paper
x_train[which(y_train == "casualty")[268], ]


## -----------------------------------------------------------------------------
classmap(vcrtrain, "survived", classCols = cols)
# classmap(vcrtrain, "survived", classCols = cols, identify = TRUE)

# red point with highest farness among highest PAC: case "d" in the paper
x_train[which(y_train == "survived")[c(14)], ] 

# near-coinciding points with highest farness: cases "e" and "f" in the paper
#
x_train[which(y_train == "survived")[265], ]

x_train[which(y_train == "survived")[287], ] 

# man --> predicted as not survived. also paid the highest 
# fare in the data and has same ticket number as passenger 
# 265. Also embarked at the same place. It turns out that
# Gustave Lesueur was the man servant of the banker Thomas Cardeza:
  
# https://www.encyclopedia-titanica.org/titanic-survivor/thomas-cardeza.html
# https://www.encyclopedia-titanica.org/titanic-survivor/gustave-lesueur.html

# blue point bottom right: case "g" in the paper
x_train[which(y_train == "survived")[90], ]

# # woman  + first class so predicted in survivors.
# # paid highest fare in whole dataset 

## -----------------------------------------------------------------------------
testdata <- data_titanic[which(data_titanic$dataType == "test"), -13]

dim(testdata)
x_test <- testdata[, -12]
y_test <- testdata[, 12]
table(y_test)

## -----------------------------------------------------------------------------
vcrtest <- vcr.rpart.newdata(x_test, y_test, vcrtrain)

confmat.vcr(vcrtest)
cols <- c("firebrick", "blue")

## -----------------------------------------------------------------------------
stackedplot(vcrtest, classCols = cols,
            main = "Stacked plot of Titanic test data")

## -----------------------------------------------------------------------------
silplot(vcrtest, classCols = cols, 
        main = "Silhouettes of rpart on Titanic test data")

## -----------------------------------------------------------------------------
classmap(vcrtest, "casualty", classCols = cols) 

## -----------------------------------------------------------------------------
classmap(vcrtest, "survived", classCols = cols) 

