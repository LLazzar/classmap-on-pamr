\name{vcr.neural.train}
\alias{vcr.neural.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare for visualization of a neural network classification on training data.
}
\description{
Produces output for the purpose of constructing graphical displays such as the \code{\link{classmap}}. The user first needs train a neural network. The representation of the data in a given layer (e.g. the final layer before applying the softmax function) then serves as the argument \code{X} to  \code{\link{vcr.neural.train}}.
}

\usage{
vcr.neural.train(X, y, probs, estmethod = meancov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the coordinates of the \code{n} objects of the training
                data, in the layer chosen by the user. Missing
                values are not allowed.
}
\item{y}{factor with the given class labels of the objects.
                Make sure that the levels are in the same order as
                used in the neural net, i.e. the columns of its
                binary "once-hot-encoded" response vectors.}
\item{probs}{posterior probabilities obtained by the neural
                net, e.g. in keras. For each case (row of \code{X}),
                the classes have probabilities that add up to 1.
                Each row of the matrix probs contains these
                probabilities. The columns of probs must be in
                the same order as the levels of \code{y}.}
\item{estmethod}{function for location and covariance estimation.
                Should return a list with \code{$m} and \code{$S}.
                Can be \code{meancov}
                (classical mean and covariance matrix) or \code{DetMCD}.
                If one or more classes have a singular covariance
                matrix, the function automatically switches to
                the PCA-based farness used in \code{\link{vcr.svm.train}}.}

}


\value{
A list with components: \cr
 \item{X}{the coordinates of the \code{n} objects of the training
                data, in the layer chosen by the user.
}
  \item{yint}{number of the given class of each case. Can contain \code{NA}'s.}
  \item{y}{given class label of each case. Can contain \code{NA}'s.}
  \item{levels}{levels of \code{y}}
  \item{predint}{predicted class number of each case. For each case this is the class with the highest posterior probability. Always exists.}
  \item{pred}{predicted label of each case.}
  \item{altint}{number of the alternative class. Among the classes different from the given class, it is the one with the highest posterior probability. Is \code{NA} for cases whose \code{y} is missing.}
  \item{altlab}{label of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
  \item{ncolX}{number of columns in \code{X}.  Keep??}
   \item{PAC}{probability of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
  \item{computeMD}{Whether or not the farness is computed using the Mahalanobis distance.}
  \item{classMS}{list with center and covariance matrix of each class}
   \item{PCAfits}{if not \code{NULL}, PCA fits to each class, estimated from the training data but also useful for new data.}
  \item{figparams}{parameters for computing \code{fig}, can be used for new data.}
  \item{fig}{distance of each case \eqn{i} from each class \eqn{g}. Always exists.}
  \item{farness}{farness of each case from its given class. Is \code{NA} for cases whose \code{y} is missing.}
  \item{ofarness}{for each case \code{i}, its lowest  \code{fig[i,g]} to any class \code{g}. Always exists.}
}
\references{
Raymaekers J., Rousseeuw P.J.(2021). Silhouettes and quasi residual plots for neural nets and tree-based classifiers. \href{https://arxiv.org/abs/2106.08814}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.neural.newdata}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}
}

\examples{
# For examples, we refer to the vignette:
\dontrun{
vignette("Neural_net_examples")
}
}
