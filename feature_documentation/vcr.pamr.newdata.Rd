\name{vcr.pamr.newdata}
\alias{vcr.pamr.newdata}
\title{
Prepare for visualization of a pamr (NSC) classification on new/test data.
}
\description{
Produces output for the purpose of constructing graphical displays such as the \code{\link{silplot}, \code{\link{qresplot} and \code{\link{mdscolorscale}}. Requires the output of \code{\link{vcr.pamr.train}} as an argument.}

\usage{
vcr.pamr.newdata(newdata, vcr.pamr.train.out, prior=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{newdata}{A list with components: \code{$x}: with genes expression levels in the rows (samples in the columns) for the new cases; \code{$y} a factor with the class labels for the new cases. The structure of \code{newdata} should mirror that of \code{data} input into \code{vcr.pamr.train.out}. It's advisable to factorize \code{newdata$y} labels in tandem with training data during initial preprocessing to avoid possible inconsistencies in the attribuition of the levels. \code{newdata}$y may include fewer unique labels than \code{data$y}.}
\item{vcr.pamr.train.out}{The output of \code{\link{vcr.pamr.train}} on the training data with \code{pamrfitcv} left to \code{NULL} (i.e. a simple evaluation on training set and not on cv results)}
\item{prior}{A different set of priors can be specifified as in \code{pamr::pamr.predict}.}
}

\value{A list with components:
  \item{yintnew}{number of the given class of each case.}
  \item{ynew}{given class label of each case.}
  \item{levels}{levels of the response, from \code{vcr.pamr.train.out}.}
  \item{predint}{predicted class number of each case. According to maximum a posteriori rule. Always exists.}
  \item{pred}{predicted label of each case.}
  \item{altint}{number of the alternative class. Among the classes different from the given class, it is the one with the highest posterior probability. Is \code{NA} for cases whose \code{y} is missing.}
  \item{altlab}{label of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
  \item{fig}{farness of each case \eqn{i} from each class \eqn{g}. Is \code{NULL} if distToClass is not feeded.}
  \item{farness}{farness of each case from its given/true class. Is \code{NULL} if distToClass is not feeded.}
  \item{ofarness}{for each case \eqn{i}, its lowest  \eqn{fig[i,g]} to any class \eqn{g}. Is \code{NULL} if distToClass is not feeded.}
  \item{threshold}{The effective threshold value evaluated.}
  \item{predictparams}{A list with the inputs used to call the \code{pamr::pamr.predict} inside \code{\link{vcr.pamr.newdata}} to get the posterior probability for the considered cases.}
}

\references{--add paper at release}


\author{
 Lazzari. L, Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.pamr.newdata}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}, \code{\link{mdscolorscale}}, \code{\link{qresplot}}
}

\examples{
# For examples, we refer to the vignette:
\dontrun{
vignette("MAKE EXAMPLE VIGNETTE")
}
}
