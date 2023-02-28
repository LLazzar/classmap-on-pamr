\name{vcr.knn.newdata}
\alias{vcr.knn.newdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Carry out a k-nearest neighbor classification on new data, and prepare to visualize its results.
}
\description{
Predicts class labels for new data by k nearest neighbors, using the output of \code{\link{vcr.knn.train}} on the training data. For cases in the new data whose given label \code{ynew} is not \code{NA}, additional output is produced for constructing graphical displays such as the \code{\link{classmap}}.
}
\usage{
vcr.knn.newdata(Xnew, ynew = NULL, vcr.knn.train.out, LOO = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xnew}{If the training data was a matrix of coordinates, \code{Xnew} must be such a matrix with the same number of columns. If the training data was a set of dissimilarities, \code{Xnew} must be a rectangular matrix of dissimilarities, with each row containing the dissmilarities of a new case to all training cases. Missing values are not allowed.}
\item{ynew}{factor with class membership of each new case. Can be \code{NA} for some or all cases. If \code{NULL}, is assumed to be \code{NA} everywhere.}
\item{vcr.knn.train.out}{output of \code{\link{vcr.knn.train}} on the training data.}
\item{LOO}{leave one out. Only used when testing this function on a subset of the training data. Default is \code{LOO=FALSE}.}
}

\value{
A list with components: \cr

\item{yintnew}{number of the given class of each case. Can contain \code{NA}'s.}
\item{ynew}{given class label of each case. Can contain \code{NA}'s.}
\item{levels}{levels of the response, from \code{vcr.knn.train.out}.}
\item{predint}{predicted class number of each case. Always exists.}
\item{pred}{predicted label of each case.}
\item{altint}{number of the alternative class. Among the classes different from the given class, it is the one with the highest posterior probability. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{altlab}{label of the alternative class. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{PAC}{probability of the alternative class. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{fig}{distance of each case \eqn{i} from each class \eqn{g}. Always exists.}
\item{farness}{farness of each case from its given class. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{ofarness}{for each case \code{i}, its lowest  \code{fig[i,g]} to any class \code{g}. Always exists.}
\item{k}{the requested number of nearest neighbors, from \code{vcr.knn.train.out}.}
\item{ktrues}{for each case this contains the actual number of elements in its neighborhood. This can be higher than \code{k} due to ties.}
\item{counts}{a matrix with 3 columns, each row representing a case. For the neighborhood of each case it says how many members it has from the given class, the predicted class, and the alternative class. The first and third entry is \code{NA} for cases whose \code{ynew} is missing.}

}

\references{
Raymaekers J., Rousseeuw P.J., Hubert M. (2021). Class maps for visualizing classification results. \emph{Technometrics}, appeared online. \doi{10.1080/00401706.2021.1927849}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.knn.train}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}
}

\examples{
data("data_floralbuds")
X <- data_floralbuds[, 1:6]; y <- data_floralbuds[, 7]
set.seed(12345); trainset <- sample(1:550, 275)
vcr.train <- vcr.knn.train(X[trainset, ], y[trainset], k = 5)
vcr.test <- vcr.knn.newdata(X[-trainset, ], y[-trainset], vcr.train)
confmat.vcr(vcr.train) # for comparison
confmat.vcr(vcr.test)
cols <- c("saddlebrown", "orange", "olivedrab4", "royalblue3")
stackedplot(vcr.train, classCols = cols) # for comparison
stackedplot(vcr.test, classCols = cols)
classmap(vcr.train, "bud", classCols = cols) # for comparison
classmap(vcr.test, "bud", classCols = cols)

# For more examples, we refer to the vignette:
\dontrun{
vignette("K_nearest_neighbors_examples")
}
}
