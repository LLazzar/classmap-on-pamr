\name{vcr.knn.train}
\alias{vcr.knn.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Carry out a k-nearest neighbor classification on training data, and prepare to visualize its results.
}
\description{
Carries out a k-nearest neighbor classification on the training data. Various additional output is produced for the purpose of constructing graphical displays such as the \code{\link{classmap}}.
}
\usage{
vcr.knn.train(X, y, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{This can be a rectangular matrix or data frame of (already standardized) measurements, or a dist object obtained from \code{\link[stats:dist]{stats::dist}} or \code{\link[cluster:daisy]{cluster::daisy}}. Missing values are not allowed.
}
\item{y}{factor with the given (observed) class labels. There need to be non-missing \code{y} in order to be able to train the classifier.}
\item{k}{the number of nearest neighbors used. It can be selected by running cross-validation using a different package.}
}

\value{
A list with components: \cr

\item{yint}{number of the given class of each case. Can contain \code{NA}'s.}
\item{y}{given class label of each case. Can contain \code{NA}'s.}
\item{levels}{levels of \code{y}}
\item{predint}{predicted class number of each case. Always exists.}
\item{pred}{predicted label of each case.}
\item{altint}{number of the alternative class. Among the classes different from the given class, it is the one with the highest posterior probability. Is \code{NA} for cases whose \code{y} is missing.}
\item{altlab}{label of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
\item{PAC}{probability of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
\item{figparams}{parameters used to compute \code{fig}.}
\item{fig}{distance of each case \eqn{i} from each class \eqn{g}. Always exists.}
\item{farness}{farness of each case from its given class. Is \code{NA} for cases whose \code{y} is missing.}
\item{ofarness}{for each case \eqn{i}, its lowest  \code{fig[i,g]} to any class \eqn{g}. Always exists.}
\item{k}{the requested number of nearest neighbors, from the arguments. Will also be used for  classifying new data.}
\item{ktrues}{for each case this contains the actual number of elements in its neighborhood. This can be higher than \code{k} due to ties.}
\item{counts}{a matrix with 3 columns, each row representing a case. For the neighborhood of each case it says how many members it has from the given class, the predicted class, and the alternative class. The first and third entry is \code{NA} for cases whose \code{y} is missing.}
\item{X}{If the argument \code{X} was a data frame or matrix of coordinates, \code{as.matrix(X)} is returned here. This is useful for classifying new data.}
}

\references{
Raymaekers J., Rousseeuw P.J., Hubert M. (2021). Class maps for visualizing classification results. \emph{Technometrics}, appeared online. \doi{10.1080/00401706.2021.1927849}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.knn.newdata}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}
}

\examples{
vcrout <- vcr.knn.train(iris[, 1:4], iris[, 5], k = 5)
confmat.vcr(vcrout)
stackedplot(vcrout)
classmap(vcrout, "versicolor", classCols = 2:4)
# The cases misclassified as virginica are shown in blue.

# For more examples, we refer to the vignette:
\dontrun{
vignette("K_nearest_neighbors_examples")
}
}
