\name{vcr.svm.train}
\alias{vcr.svm.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare for visualization of a support vector machine classification on training data.
}
\description{
Produces output for the purpose of constructing graphical displays such as the \code{\link{classmap}}. The user first needs to run a support vector machine classification on the data by \code{\link[e1071:svm]{e1071::svm}}, with the option \code{probability = TRUE}. This classification can be with two or more classes. The output of \code{\link[e1071:svm]{e1071::svm}} is then an argument to \code{\link{vcr.svm.train}}. As \code{\link[e1071:svm]{e1071::svm}} does not output the data itself, it needs to be given as well, in the arguments \code{X} and \code{y}.
}
\usage{
vcr.svm.train(X, y, svfit, ortho = FALSE)
}

\arguments{

\item{X}{matrix of data coordinates, as used in \code{\link[e1071:svm]{e1071::svm}}. Missing values are not allowed.}
\item{y}{factor with the given (observed) class labels. It is crucial that X and y are exactly the same as in the call to \code{\link[e1071:svm]{e1071::svm}}.}
\item{svfit}{an object returned by \code{\link[e1071:svm]{e1071::svm}}, called with exactly the same \code{X} and \code{y} as above.}
\item{ortho}{If \code{TRUE}, will compute farness in the orthogonal complement of the vector beta given by \code{\link[e1071:svm]{e1071::svm}}. Is only possible for 2 classes, else there would be several beta vectors.}
}

\value{
A list with components: \cr
\item{yint}{number of the given class of each case. Can contain \code{NA}'s.}
\item{y}{given class label of each case. Can contain \code{NA}'s.}
\item{levels}{levels of the response \code{y}.}
\item{predint}{predicted class number of each case. Always exists.}
\item{pred}{predicted label of each case.}
\item{altint}{number of the alternative class. Among the classes different from the given class, it is the one with the highest posterior probability. Is \code{NA} for cases whose \code{y} is missing.}
\item{altlab}{label of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
\item{PAC}{probability of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
\item{figparams}{parameters used in \code{fig}, can be used for new data.}
\item{fig}{distance of each case \eqn{i} from each class \eqn{g}. Always exists.}
\item{farness}{farness of each case from its given class. Is \code{NA} for cases whose \code{y} is missing.}
\item{ofarness}{for each case \eqn{i}, its lowest \code{fig[i,g]} to any class \code{g}. Always exists.}
\item{svfit}{as it was input, will be useful for new data.}
\item{X}{the matrix of data coordinates from the arguments. This is useful for classifying new data.}
}

\references{
Raymaekers J., Rousseeuw P.J., Hubert M. (2021). Class maps for visualizing classification results. \emph{Technometrics}, appeared online. \doi{10.1080/00401706.2021.1927849}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.knn.newdata}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}, \code{\link[e1071:svm]{e1071::svm}}
}

\examples{
library(e1071)
set.seed(1); X <- matrix(rnorm(200 * 2), ncol = 2)
X[1:100, ] <- X[1:100, ] + 2
X[101:150, ] <- X[101:150, ] - 2
y <- as.factor(c(rep("blue", 150), rep("red", 50)))
cols <- c("deepskyblue3", "red")
plot(X, col = cols[as.numeric(y)], pch = 19)
# We now fit an SVM with radial basis kernel to the data:
set.seed(1) # to make the result of svm() reproducible.
svmfit <- svm(y~., data = data.frame(X = X, y = y),
scale = FALSE, kernel = "radial", cost = 10,
gamma = 1, probability = TRUE)
plot(svmfit$decision.values, col = cols[as.numeric(y)]); abline(h = 0)
# so the decision values separate the classes reasonably well.
plot(svmfit, data = data.frame(X = X, y = y), X.2~X.1, col = cols)
# The boundary is far from linear (but in feature space it is).
vcr.train <- vcr.svm.train(X, y, svfit = svmfit)
confmat.vcr(vcr.train)
stackedplot(vcr.train, classCols = cols)
classmap(vcr.train, "blue", classCols = cols)
classmap(vcr.train, "red", classCols = cols)

# For more examples, we refer to the vignette:
\dontrun{
vignette("Support_vector_machine_examples")
}
}
