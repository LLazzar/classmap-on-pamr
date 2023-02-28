\name{vcr.rpart.newdata}
\alias{vcr.rpart.newdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare for visualization of an rpart classification on new data.
}
\description{
Produces output for the purpose of constructing graphical displays such as the \code{\link{classmap}} on new data. Requires the output of
\code{\link{vcr.rpart.train}} as an argument.
}

\usage{
vcr.rpart.newdata(Xnew, ynew = NULL, vcr.rpart.train.out,
                  LOO = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xnew}{data matrix of the new data, with the same
  number of columns \code{d} as in the training data. Missing values are not allowed.
}
\item{ynew}{factor with class membership of each new case. Can be \code{NA} for some or all cases. If \code{NULL}, is assumed to be \code{NA} everywhere.}
\item{vcr.rpart.train.out}{output of \code{\link{vcr.rpart.train}} on the training data.}
\item{LOO}{leave one out. Only used when testing this function on a subset of the training data. Default is \code{LOO=FALSE}.}
}


\value{
A list with components: \cr
 \item{yintnew}{number of the given class of each case. Can contain \code{NA}'s.}
\item{ynew}{given class label of each case. Can contain \code{NA}'s.}
\item{levels}{levels of the response, from \code{vcr.rpart.train.out}.}
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
\code{\link{vcr.rpart.train}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}
}

\examples{
library(rpart)
data("data_titanic")
traindata <- data_titanic[which(data_titanic$dataType == "train"), -13]
str(traindata); table(traindata$y)
set.seed(123) # rpart is not deterministic
rpart.out <- rpart(y ~ Pclass + Sex + SibSp +
                    Parch + Fare + Embarked,
                  data = traindata, method = 'class', model = TRUE)
y_train <- traindata[, 12]
x_train <- traindata[, -12]
mytype <- list(nominal = c("Name", "Sex", "Ticket", "Cabin", "Embarked"), ordratio = c("Pclass"))
# These are 5 nominal columns, and one ordinal.
# The variables not listed are by default interval-scaled.
vcrtrain <- vcr.rpart.train(x_train, y_train, rpart.out, mytype)
testdata <- data_titanic[which(data_titanic$dataType == "test"), -13]
dim(testdata)
x_test <- testdata[, -12]
y_test <- testdata[, 12]
vcrtest <- vcr.rpart.newdata(x_test, y_test, vcrtrain)
confmat.vcr(vcrtest)
silplot(vcrtest, classCols = c(2, 4))
classmap(vcrtest, "casualty", classCols = c(2, 4))
classmap(vcrtest, "survived", classCols = c(2, 4))

# For more examples, we refer to the vignette:
\dontrun{
vignette("Rpart_examples")
}
}
