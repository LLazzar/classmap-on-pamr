\name{vcr.da.newdata}
\alias{vcr.da.newdata}

\title{
Carry out discriminant analysis on new data, and prepare to visualize its results.
}
\description{
Predicts class labels for new data by discriminant analysis, using the output of \code{\link{vcr.da.train}} on the training data. For new data cases whose label in \code{yintnew} is non-missing, additional output is produced for constructing graphical displays such as the \code{\link{classmap}}.
}
\usage{
vcr.da.newdata(Xnew, ynew=NULL, vcr.da.train.out)
}

\arguments{
  \item{Xnew}{data matrix of the new data, with the same number of columns as in the training data. Missing values are not allowed.
}
\item{ynew}{factor with class membership of each new case. Can be \code{NA} for some or all cases. If \code{NULL}, is assumed to be \code{NA} everywhere.}
\item{vcr.da.train.out}{output of \code{\link{vcr.da.train}} on the training data.}
}

\value{
A list with components: \cr

\item{yintnew}{number of the given class of each case. Can contain \code{NA}'s.}
\item{ynew }{given class label of each case. Can contain \code{NA}'s.}
\item{levels }{levels of the response, from \code{vcr.da.train.out}.}
\item{predint}{predicted class number of each case. Always exists.}
\item{pred}{predicted label of each case.}
\item{altint}{number of the alternative class. Among the classes different from the given class, it is the one with the highest posterior probability. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{altlab}{label of the alternative class. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{PAC}{probability of the alternative class. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{fig}{distance of each case \eqn{i} to each class \eqn{g}. Always exists.}
\item{farness}{farness of each case \eqn{i} from its given class. Is \code{NA} for cases whose \code{ynew} is missing.}
\item{ofarness}{For each case \eqn{i}, its lowest \code{fig[i,g]} to any class \code{g}. Always exists.}
\item{classMS}{list with center and covariance matrix of each class, from \code{vcr.da.train.out}.}
\item{lCurrent}{log of mixture density of each case in its given class. Is \code{NA} for cases with missing \code{ynew}.}
\item{lPred}{log of mixture density of each case in its predicted class. Always exists.}
\item{lAlt}{log of mixture density of each case in its alternative class. Is \code{NA} for cases with missing \code{ynew}.}
}
\references{
Raymaekers J., Rousseeuw P.J., Hubert M. (2021). Class maps for visualizing classification results. \emph{Technometrics}, appeared online. \doi{10.1080/00401706.2021.1927849}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.da.train}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}
}

\examples{

vcr.train <- vcr.da.train(iris[, 1:4], iris[, 5])
inds <- c(51:150) # a subset, containing only 2 classes
iris2 <- iris[inds, ] # fake "new" data
iris2[c(1:10, 51:60), 5] <- NA
vcr.test <- vcr.da.newdata(iris2[, 1:4], iris2[, 5], vcr.train)
vcr.test$PAC[1:25] # between 0 and 1. Is NA where the response is.
plot(vcr.test$PAC, vcr.train$PAC[inds]); abline(0, 1) # match
plot(vcr.test$farness, vcr.train$farness[inds]); abline(0, 1) # match
confmat.vcr(vcr.train) # for comparison
confmat.vcr(vcr.test)
stackedplot(vcr.train) # for comparison
stackedplot(vcr.test)
classmap(vcr.train, "versicolor", classCols = 2:4) # for comparison
classmap(vcr.test, "versicolor", classCols = 2:4) # has fewer points

# For more examples, we refer to the vignette:
\dontrun{
vignette("Discriminant_analysis_examples")
}
}

