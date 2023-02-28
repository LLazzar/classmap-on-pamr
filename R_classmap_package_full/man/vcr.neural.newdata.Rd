\name{vcr.neural.newdata}
\alias{vcr.neural.newdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare for visualization of a neural network classification on new data.
}
\description{
  Prepares graphical display of new data fitted by a neural
  net that was modeled on the training data, using the output
  of \code{\link{vcr.neural.train}} on the training data.
}
\usage{
vcr.neural.newdata(Xnew, ynew = NULL, probs,
                   vcr.neural.train.out)
}

\arguments{
\item{Xnew}{data matrix of the new data, with the same number of columns as in the training data. Missing values in \code{Xnew} are not allowed.}
\item{ynew}{factor with class membership of each new case. Can be \code{NA} for some or all cases. If \code{NULL}, is assumed to be \code{NA} everywhere.}
\item{probs}{posterior probabilities obtained by running the neural net on the new data.}
\item{vcr.neural.train.out}{ output of \code{\link{vcr.neural.train}} on the training data.}
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
Raymaekers J., Rousseeuw P.J.(2021). Silhouettes and quasi residual plots for neural nets and tree-based classifiers. \href{https://arxiv.org/abs/2106.08814}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.neural.train}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}
}

\examples{
# For examples, we refer to the vignette:
\dontrun{
vignette("Neural_net_examples")
}
}
