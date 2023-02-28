\name{vcr.forest.newdata}
\alias{vcr.forest.newdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare for visualization of a random forest classification on new data.
}
\description{
Produces output for the purpose of constructing graphical displays such as the \code{\link{classmap}} on new data. Requires the output of
\code{\link{vcr.forest.train}} as an argument.
}

\usage{
vcr.forest.newdata(Xnew, ynew = NULL, vcr.forest.train.out,
                   LOO = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xnew}{data matrix of the new data, with the same
                         number of columns \code{d} as in the training data.
                         Missing values are not allowed.
}
\item{ynew}{factor with class membership of each new case. Can be \code{NA} for some or all cases. If \code{NULL}, is assumed to be \code{NA} everywhere.}
\item{vcr.forest.train.out}{output of \code{\link{vcr.forest.train}} on the training data.}
\item{LOO}{leave one out. Only used when testing this function on a subset of the training data. Default is \code{LOO=FALSE}.}
}


\value{
A list with components: \cr
 \item{yintnew}{number of the given class of each case. Can contain \code{NA}'s.}
\item{ynew}{given class label of each case. Can contain \code{NA}'s.}
\item{levels}{levels of the response, from \code{vcr.forest.train.out}.}
\item{predint}{predicted class number of each case. Always exists.}
\item{pred}{predicted label of each case.}
\item{altint}{number of the alternative class. Among the classes different from the given class, it is the one with the highest posterior probability. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{altlab}{alternative label if yintnew was given, else \code{NA}.}
\item{PAC}{probability of the alternative class. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{fig}{distance of each case \eqn{i} from each class \eqn{g}. Always exists.}
\item{farness}{farness of each case from its given class. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{ofarness}{for each case \eqn{i}, its lowest \code{fig[i,g]} to any class \eqn{g}. Always exists.}
}
\references{
Raymaekers J., Rousseeuw P.J.(2021). Silhouettes and quasi residual plots for neural nets and tree-based classifiers. \href{https://arxiv.org/abs/2106.08814}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.forest.train}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}
}

\examples{
library(randomForest)
data("data_instagram")
traindata <- data_instagram[which(data_instagram$dataType == "train"), -13]
set.seed(71) # randomForest is not deterministic
rfout <- randomForest(y ~ ., data = traindata, keep.forest = TRUE)
mytype <- list(symm = c(1, 5, 7, 8)) # These 4 columns are
# (symmetric) binary variables. The variables that are not
# listed are interval-scaled by default.
x_train <- traindata[, -12]
y_train <- traindata[, 12]
vcrtrain <- vcr.forest.train(X = x_train, y = y_train,
                            trainfit = rfout, type = mytype)
testdata <- data_instagram[which(data_instagram$dataType == "test"), -13]
Xnew <- testdata[, -12]
ynew <- testdata[, 12]
vcrtest <- vcr.forest.newdata(Xnew, ynew, vcrtrain)
confmat.vcr(vcrtest)
stackedplot(vcrtest, classCol = c(4, 2))
silplot(vcrtest, classCols = c(4, 2))
classmap(vcrtest, "genuine", classCols = c(4, 2))
classmap(vcrtest, "fake", classCols = c(4, 2))

# For more examples, we refer to the vignette:
\dontrun{
vignette("Random_forest_examples")
}
}
