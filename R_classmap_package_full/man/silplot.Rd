\name{silplot}
\alias{silplot}
\title{
Draw the silhouette plot of a classification
}
\description{
Draw the silhouette plot to visualize classification results, based on the output of one of the \code{vcr.*.*} functions in this package. The horizontal axis of the silhouette plot shows each case's \code{s(i)}.
}
\usage{
silplot(vcrout, classLabels = NULL, classCols = NULL,
        showLegend = TRUE, showClassNumbers = FALSE,
        showCases = FALSE, drawLineAtAverage = FALSE,
        topdown = TRUE, main = NULL, summary = TRUE)
}

\arguments{
 \item{vcrout}{output of \code{vcr.*.train} or \code{vcr.*.newdata}. Required.}
 \item{classLabels}{the labels (levels) of the classes. If \code{NULL}, they are taken from \code{vcrout}.}
 \item{classCols}{a list of colors for the classes. There should be at least as many as there are levels. If \code{NULL} a default palette is used.}

\item{showLegend}{if \code{TRUE}, a legend is shown to the right of the plot.}
\item{showClassNumbers}{if \code{TRUE}, the legend will show the class numbers
 instead of the class labels.}
\item{showCases}{if \code{TRUE}, the plot shows the numbers of the cases.
 They are only readable when the number of cases
 is relatively small.}
\item{topdown}{if \code{TRUE} (the default), the silhouettes are
 plotted from top to bottom. Otherwise they
 are plotted from left to right.}
\item{drawLineAtAverage}{if \code{TRUE}, drwas a line at the average value
 of the \code{s(i)}.}
\item{main }{title for the plot. If \code{NULL}, a default title is used.}
\item{summary}{if \code{TRUE}, puts a summary table on the screen
 with for each class its number, label, number of class members,
 and the average of its \code{s(i)}.}

}

\value{
A ggplot object containing the silhouette plot.
}

\references{
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
silplot(vcrout)
# For more examples, we refer to the vignettes:
\dontrun{
vignette("Discriminant_analysis_examples")
vignette("K_nearest_neighbors_examples")
vignette("Support_vector_machine_examples")
vignette("Rpart_examples")
vignette("Forest_examples")
vignette("Neural_net_examples")
}
}

