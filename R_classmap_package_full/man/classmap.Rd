\name{classmap}
\alias{classmap}
\title{
Draw the class map to visualize classification results.
}
\description{
Draw the class map to visualize classification results, based on the output of one of the
\code{vcr.*.*} functions in this package. The vertical axis of the class map shows each case's \code{PAC}, the conditional probability that it belongs to an alternative class. The \code{farness} on the horizontal axis is the probability of a member of the given class being at most as far from the class as the case itself.
}
\usage{
classmap(vcrout, whichclass, classLabels = NULL, classCols = NULL,
         main = NULL, cutoff = 0.99, plotcutoff = TRUE,
         identify = FALSE, cex = 1, cex.main = 1.2, cex.lab = NULL,
         cex.axis = NULL, opacity = 1,
         squareplot = TRUE,  maxprob = NULL, maxfactor = NULL)
}

\arguments{
 \item{vcrout}{output of \code{vcr.*.train} or \code{vcr.*.newdata}. Required.}

\item{whichclass}{the number or level of the class to be displayed. Required.}
 \item{classLabels}{the labels (levels) of the classes. If \code{NULL}, they are taken from \code{vcrout}.}
 \item{classCols}{a list of colors for the class labels. There should be at least as many as there are levels. If \code{NULL} the \code{classCols} are taken as 2, 3, 4, ...}
 \item{main}{title for the plot.}
  \item{cutoff}{cases with overall farness \code{vcrout$ofarness} > \code{cutoff} are flagged as outliers.}
   \item{plotcutoff}{If true, plots the cutoff on the farness values as a vertical line.}
 \item{identify}{if \code{TRUE}, left-click on a point to get its number, then ESC to exit.}
 \item{cex}{passed on to \code{\link[graphics:plot]{graphics::plot}}.}
 \item{cex.main}{same, for title.}
 \item{cex.lab}{same, for labels on horizontal and vertical axes.}
 \item{cex.axis}{same, for axes.}
 \item{opacity}{determines opacity of plotted dots. Value between 0 and 1, where 0 is transparent and 1 is opaque.}
 \item{squareplot}{If \code{TRUE}, makes the axes of the plot equally long.}
 \item{maxprob}{draws the farness axis at least upto probability maxprob. If \code{NULL}, the limits are obtained automatically.}
 \item{maxfactor}{if not \code{NULL}, a number slightly higher than 1 to increase the space at the right hand side of the plot, to make room for marking points.}
}

\value{
Executing the function plots the class map and returns
\item{coordinates}{ a matrix with 2 columns containing the coordinates of the plotted points. The first coordinate is the quantile of the farness probability. This makes it easier to add text next to interesting points. If \code{identify = T}, the attribute \code{ids} of \code{coordinates} contains the row numbers of the identified points in the matrix \code{coordinates}.}
}

\references{
Raymaekers J., Rousseeuw P.J., Hubert M. (2021). Class maps for visualizing classification results. \emph{Technometrics}, appeared online. \doi{10.1080/00401706.2021.1927849}{(link to open access pdf)}

Raymaekers J., Rousseeuw P.J.(2021). Silhouettes and quasi residual plots for neural nets and tree-based classifiers. \href{https://arxiv.org/abs/2106.08814}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.da.train}}, \code{\link{vcr.da.newdata}},\cr \code{\link{vcr.knn.train}}, \code{\link{vcr.knn.newdata}},\cr \code{\link{vcr.svm.train}}, \code{\link{vcr.svm.newdata}},\cr \code{\link{vcr.rpart.train}}, \code{\link{vcr.rpart.newdata}},\cr \code{\link{vcr.forest.train}}, \code{\link{vcr.forest.newdata}},\cr \code{\link{vcr.neural.train}}, \code{\link{vcr.neural.newdata}}
}

\examples{
vcrout <- vcr.da.train(iris[, 1:4], iris[, 5])
classmap(vcrout, "setosa", classCols = 2:4) # tight class
classmap(vcrout, "versicolor", classCols = 2:4) # less tight
# The cases misclassified as virginica are shown in blue.
classmap(vcrout, "virginica", classCols = 2:4)
# The case misclassified as versicolor is shown in green.

# For more examples, we refer to the vignettes:
\dontrun{
vignette("Discriminant_analysis_examples")
vignette("K_nearest_neighbors_examples")
vignette("Support_vector_machine_examples")
vignette("Rpart_examples")
vignette("Random_forest_examples")
vignette("Neural_net_examples")
}
}

