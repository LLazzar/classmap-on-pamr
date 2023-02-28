\name{vcr.da.train}
\alias{vcr.da.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Carry out discriminant analysis on training data, and prepare to visualize its results.
}
\description{
Custom DA function which prepares for graphical displays such as the \code{\link{classmap}}. The disciminant analysis itself is carried out by the maximum a posteriori rule, which maximizes the density of the mixture.
}
\usage{
vcr.da.train(X, y, rule = "QDA", estmethod = "meancov")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numerical matrix containing the predictors in its columns. Missing values are not allowed.
}
\item{y}{a factor with the given class labels.}
\item{rule}{either "\code{QDA}" for quadratic discriminant analysis or "\code{LDA}"  for linear discriminant analysis.}
\item{estmethod}{function for location and covariance estimation.
Should return a list with the center \code{$m} and the covariance matrix \code{$S}. The default is \code{"meancov"} (classical mean and covariance matrix), and the option \code{"DetMCD"} (based on \code{\link[robustbase:covMcd]{robustbase::covMcd}}) is also provided.}
}


\value{
A list with components: \cr
  \item{yint}{number of the given class of each case. Can contain \code{NA}'s.}
  \item{y}{given class label of each case. Can contain \code{NA}'s.}
  \item{levels}{levels of \code{y}}
  \item{predint}{predicted class number of each case. For each case this is the class with the highest posterior probability. Always exists.}
  \item{pred}{predicted label of each case.}
  \item{altint}{number of the alternative class. Among the classes different from the given class, it is the one with the highest posterior probability. Is \code{NA} for cases whose \code{y} is missing.}
  \item{altlab}{label of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
  \item{PAC}{probability of the alternative class. Is \code{NA} for cases whose \code{y} is missing.}
  \item{figparams}{parameters for computing \code{fig}, can be used for new data.}
  \item{fig}{distance of each case \eqn{i} from each class \eqn{g}. Always exists.}
  \item{farness}{farness of each case from its given class. Is \code{NA} for cases whose \code{y} is missing.}
  \item{ofarness}{for each case \code{i}, its lowest  \code{fig[i,g]} to any class \code{g}. Always exists.}
  \item{classMS}{list with center and covariance matrix of each class}
  \item{lCurrent}{log of mixture density of each case in its given class. Is \code{NA} for cases with missing \code{y}.
}
  \item{lPred}{log of mixture density of each case in its predicted class. Always exists.
}
  \item{lAlt}{log of mixture density of each case in its alternative class. Is \code{NA} for cases with missing \code{y}.}
}
\references{
Raymaekers J., Rousseeuw P.J., Hubert M. (2021). Class maps for visualizing classification results. \emph{Technometrics}, appeared online. \doi{10.1080/00401706.2021.1927849}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.da.newdata}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}
}

\examples{

data("data_floralbuds")
X <- data_floralbuds[, 1:6]; y <- data_floralbuds[, 7]
vcrout <- vcr.da.train(X, y, rule = "QDA")
# For linear discriminant analysis, put rule = "LDA".
confmat.vcr(vcrout) # There are a few outliers
cols <- c("saddlebrown", "orange", "olivedrab4", "royalblue3")
stackedplot(vcrout, classCols = cols)
classmap(vcrout, "bud", classCols = cols)

# For more examples, we refer to the vignette:
\dontrun{
vignette("Discriminant_analysis_examples")
}
}
