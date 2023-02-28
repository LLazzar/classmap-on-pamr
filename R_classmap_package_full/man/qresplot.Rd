\name{qresplot}
\alias{qresplot}

\title{
Draw a quasi residual plot of PAC versus a data feature
}
\description{
Draw a quasi residual plot to visualize classification results.
The vertical axis of the quasi residual plot shows each case's probability of alternative class (PAC). The horizontal axis
shows the feature given as the second argument in the function call.
}
\usage{
qresplot(PAC, feat, xlab = NULL, xlim = NULL,
         main = NULL, identify = FALSE, gray = TRUE,
         opacity = 1, squareplot = FALSE, plotLoess = FALSE,
         plotErrorBars = FALSE, plotQuantiles = FALSE,
         grid = NULL, probs = c(0.5, 0.75),
         cols = NULL, fac = 1, cex = 1,
         cex.main = 1.2, cex.lab = 1,
         cex.axis = 1, pch = 19)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{PAC}{vector with the PAC values of a classification,
typically the \code{$PAC} in the return of a call to a function \code{vcr.*.*}}
\item{feat}{the PAC will be plotted versus this data
feature. Note that feat does not have to be
one of the explanatory variables of the model.
It can be another variable, a combination of
variables (like a sum or a principal component
score), the row number of the cases if they
were recorded succesively, etc.}
\item{xlab}{label for the horizontal axis, i.e. the name
of variable feat.}
\item{xlim}{limits for the horizontal axis. If \code{NULL}, the
range of feat is used.}
\item{main}{title for the plot.}
\item{identify }{if \code{TRUE}, left-click on a point to get its number, then ESC to exit.}
\item{gray}{logical, if \code{TRUE} (the default) the plot region
where \eqn{PAC < 0.5} gets a light gray background.
Points in this region were classified into
their given class, and the points above this
region were misclassified.}
\item{opacity}{determines opacity of plotted dots.
Value between 0 and 1, where 0 is transparent
and 1 is opaque.}
\item{squareplot}{if \code{TRUE}, the horizontal and vertical axis will
get the same length.}
\item{plotLoess}{if \code{TRUE}, a standard loess curve is fitted and
superimposed on the plot. May not work well if
feat is discrete with few values.
At most one of the options \code{plotLoess},
\code{plotErrorbars}, or \code{plotQuantiles} can be selected.}
\item{plotErrorBars}{if \code{TRUE}, the average PAC and its standard
error are computed on the intervals of a grid
(see option grid). Then a red curve connecting
the averages is plotted, as well as two blue
curves corresponding to the average plus or
minus one standard error. At most one of the options \code{plotLoess},
\code{plotErrorbars}, or \code{plotQuantiles} can be selected.}
\item{plotQuantiles}{if \code{TRUE}, one or more quantiles of the PAC
are computed on the intervals of a grid
(see option grid). The quantiles correspond
the probabilities in option probs.
Then the curves connecting the quantiles
are plotted. At most one of the options \code{plotLoess},
\code{plotErrorbars}, or \code{plotQuantiles} can be selected.}
\item{grid}{only used when \code{plotErrorBars} or \code{plotQuantiles}
are selected. This is a vector with increasing
feat values, forming the grid. If \code{NULL}, the
grid consists of the minimum and the maximum
of feat, with 9 equispaced points between them.}
\item{probs}{only used when \code{plotQuantiles} is selected. This
is a vector with probabilities determining the
quantiles. If \code{NULL}, defaults to \code{c(0.5, 0.75)}.}
\item{cols}{only used when plotquantiles is selected.
A vector with the colors of the quantile curves.
If \code{NULL} the cols are taken as 2, 3, ...}
\item{fac}{only used when \code{plotLoess}, \code{plotErrorBars} or
\code{plotQuantiles} are selected. A real number to
multiply the resulting curves. A value \code{fac > 1}
can be useful to better visualize the curves
when they would be too close to zero.
By default (\code{fac = 1}) this is not done.}
\item{cex}{passed on to \code{\link[graphics]{plot}}.}
\item{cex.main}{same, for title.}
\item{cex.lab}{same, for labels on horizontal and vertical axes.}
\item{cex.axis}{same, for axes.}
\item{pch}{plot character for the points, defaults to 19.}

}


\value{

\item{coordinates}{a matrix with 2 columns containing the
coordinates of the plotted points. This makes it
easier to add text next to interesting points.
If \code{identify = TRUE}, the attribute ids of coordinates
contains the row numbers of the identified points
in the matrix coordinates.

}
}
\references{
Raymaekers J., Rousseeuw P.J.(2021). Silhouettes and quasi residual plots for neural nets and tree-based classifiers. \href{https://arxiv.org/abs/2106.08814}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}


\examples{
library(rpart)
data("data_titanic")
traindata <- data_titanic[which(data_titanic$dataType == "train"), -13]
set.seed(123) # rpart is not deterministic
rpart.out <- rpart(y ~ Pclass + Sex + SibSp +
                    Parch + Fare + Embarked,
                  data = traindata, method = 'class', model = TRUE)
mytype <- list(nominal = c("Name", "Sex", "Ticket", "Cabin", "Embarked"), ordratio = c("Pclass"))
x_train <- traindata[, -12]
y_train <- traindata[,  12]
vcrtrain <- vcr.rpart.train(x_train, y_train, rpart.out, mytype)
# Quasi residual plot versus age, for males only:
PAC <- vcrtrain$PAC[which(x_train$Sex == "male")]
feat <- x_train$Age[which(x_train$Sex == "male")]
qresplot(PAC, feat, xlab = "Age (years)", opacity = 0.5,
         main = "quasi residual plot for male passengers",
         plotLoess = TRUE)
text(x = 14, y = 0.60, "loess curve", col = "red", cex = 1)
}
