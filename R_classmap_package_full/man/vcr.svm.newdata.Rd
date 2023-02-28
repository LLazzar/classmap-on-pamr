\name{vcr.svm.newdata}
\alias{vcr.svm.newdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare for visualization of a support vector machine classification on new data.
}
\description{
Carries out a support vector machine classification of new data using the output of \code{\link{vcr.svm.train}} on the training data, and computes the quantities needed for its visualization.
}
\usage{
vcr.svm.newdata(Xnew, ynew = NULL, vcr.svm.train.out)
}

\arguments{

\item{Xnew}{data matrix of the new data, with the same number of columns as in the training data. Missing values in \code{Xnew} are not allowed.}
\item{ynew}{factor with class membership of each new case. Can be \code{NA} for some or all cases. If \code{NULL}, is assumed to be \code{NA} everywhere.}
\item{vcr.svm.train.out}{ output of \code{\link{vcr.svm.train}} on the training data.
}
}

\value{
A list with components: \cr
\item{yintnew}{number of the given class of each case. Can contain \code{NA}'s.}
\item{ynew}{given class label of each case. Can contain \code{NA}'s.}
\item{levels}{levels of the response, from \code{vcr.svm.train.out}.}
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
Raymaekers J., Rousseeuw P.J., Hubert M. (2021). Class maps for visualizing classification results. \emph{Technometrics}, appeared online. \doi{10.1080/00401706.2021.1927849}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.svm.train}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}, \code{\link[e1071:svm]{e1071::svm}}
}

\examples{

library(e1071)
set.seed(1); X <- matrix(rnorm(200 * 2), ncol = 2)
X[1:100, ] <- X[1:100, ] + 2
X[101:150, ] <- X[101:150, ] - 2
y <- as.factor(c(rep("blue", 150), rep("red", 50)))
# We now fit an SVM with radial basis kernel to the data:
set.seed(1) # to make the result of svm() reproducible.
svmfit <- svm(y~., data = data.frame(X = X, y = y),
scale = FALSE, kernel = "radial", cost = 10,
gamma = 1, probability = TRUE)
vcr.train <- vcr.svm.train(X, y, svfit = svmfit)
# As "new" data we take a subset of the training data:
inds <- c(1:25, 101:125, 151:175)
vcr.test <- vcr.svm.newdata(X[inds, ], y[inds], vcr.train)
plot(vcr.test$PAC, vcr.train$PAC[inds]); abline(0, 1) # match
plot(vcr.test$farness, vcr.train$farness[inds]); abline(0, 1)
confmat.vcr(vcr.test)
cols <- c("deepskyblue3", "red")
stackedplot(vcr.test, classCols = cols)
classmap(vcr.train, "blue", classCols = cols) # for comparison
classmap(vcr.test, "blue", classCols = cols)
classmap(vcr.train, "red", classCols = cols) # for comparison
classmap(vcr.test, "red", classCols = cols)


# For more examples, we refer to the vignette:
\dontrun{
vignette("Support_vector_machine_examples")
}
}
