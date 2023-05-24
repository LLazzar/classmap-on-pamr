\name{vcr.pamr.train}
\alias{vcr.pamr.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare for visualization of a pamr (NSC) classification on training data.
}
\description{
Produces output for the purpose of constructing graphical displays such as the \code{silplot}, \code{qresplot} and \code{\link{mdscolorscale}}. The user first needs train a pamr model (Nearest Shrunken Centroid) on the data by \code{pamr::train}. This then serves as an argument to \code{\link{vcr.pamr.train}}.
}

\usage{
vcr.pamr.train(data, pamrfit, pamrfitcv=NULL, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{data}{Exactly the same input data used in pamr.train that is a list with components: $x: with genes expression levels in the rows (samples in the columns); $y is required to be factor with the class labels for each sample. Data should be feeded to \code{pamr.train} already factorized. Additional components: $genenames, a vector of gene names, and $geneid a vector of gene identifiers.}
\item{pamrfit}{The output of a \code{pamr::pamr.train} run}
\item{pamrfitcv}{The output of a \code{pamr::pamr.cv} run with the \code{pamrfit} object. If feeded the cross-validated results are prepared for visualization instead of the results on the normal training process. Default is \code{NULL}}
\item{threshold}{Define the threshold at which training classification results are to be evaluated. This is because \code{pamrfit} includes multiple models, each linked to a unique threshold from a vector of discrete values. The first threshold value in \code{pamrfit}, exceeding the inputted \code{threshold}, is selected.}
}


\value{
A list with components: \cr
  \item{yint}{number of the given class of each case.}
  \item{y}{given class label of each case.}
  \item{levels}{levels of \code{y}}
  \item{predint}{predicted class number of each case. For each case this is the class with the highest posterior probability. Always exists.}
  \item{pred}{predicted label of each case.}
  \item{altint}{number of the alternative class. Among the classes different from the given class, it is the one with the highest posterior probability. Is \code{NA} for cases whose \code{y} is missing.}
  \item{altlab}{label of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
  \item{fig}{farness of each case \eqn{i} from each class \eqn{g}. Is \code{NULL} if distToClass is not feeded.}
  \item{farness}{farness of each case from its given/true class. Is \code{NULL} if distToClass is not feeded.}
  \item{ofarness}{for each case \code{i}, its lowest  \code{fig[i,g]} to any class \code{g}. Is \code{NULL} if distToClass is not feeded.}
}
  \item{pamrfit}{The \code{pamrfit} object taken as input.}
  \item{pamrfitcv}{The \code{pamrfitcv} object taken as input, if provided.}
  \item{threshold}{The effective threshold value evaluated.}
  \item{ii}{The index of the threshold value evaluated.}
}
\references{--add paper at release}
}

\author{
 Lazzari. L, Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.pamr.newdata}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}, \code{\link{mdscolorscale}, \code{\link{qresplot}}
}

\examples{
# For examples, we refer to the vignette:
\dontrun{
vignette("MAKE EXAMPLE VIGNETTE")
}
}
