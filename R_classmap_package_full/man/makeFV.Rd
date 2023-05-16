\name{makeFV}
\alias{makeFV}
\title{
Constructs feature vectors from a kernel matrix.
}
\description{
Constructs feature vectors from a kernel matrix.
}
\usage{
makeFV(kmat, transfmat = NULL, precS = 1e-12)
}

\arguments{

\item{kmat}{a kernel matrix. If \code{transfmat} is \code{NULL}, we are dealing with training data and then \code{kmat} must be a square kernel matrix (of size \eqn{n} by \eqn{n} when there are \eqn{n} cases). Such a PSD matrix kmat can e.g. be produced by \code{\link{makeKernel}} or by \code{\link[kernlab:kernelMatrix]{kernlab::kernelMatrix}}. If on the other hand \code{transfmat} is not \code{NULL}, we are dealing with a test set. See details for the precise working.}
\item{transfmat}{transformation matrix. If not \code{NULL}, it is the value \code{transfmat} of \code{\link{makeFV}} on training data. It has to be a square matrix, with as many rows as there were training data.}
\item{precS}{if not \code{NULL}, eigenvalues of \code{kmat} below \code{precS} will be set equal to precS.}
}

\details{
If \code{transfmat} is non-\code{NULL}, we are dealing with a test set.
Denote the number of cases in the test set by \eqn{m \geq 1}. Each row of \code{kmat} of the test set then must contain the kernel values of a new case with all cases in the training set. Therefore the kernel matrix kmat must have dimensions \eqn{m} by \eqn{n}. The matrix \code{kmat} can e.g. be produced by \code{\link{makeKernel}}. It can also be obtained by running \code{\link[kernlab:kernelMatrix]{kernlab::kernelMatrix}} on the union of the training set and the test set, yielding an \eqn{(n+m)} by \eqn{(n+m)} matrix, from which one then takes the \eqn{[(n+1):m , 1:n]} submatrix.
}


\value{
A list with components: \cr
\item{Xf}{When makeKV is applied to the training set, \code{Xf} has coordinates of \eqn{n} points (vectors), the plain inner products of which equal the kernel matrix of the training set. That is, \code{kmat} = \code{Xf} \code{Xf}'. The \code{Xf} are expressed in an orthogonal basis in which the variance of the coordinates is decreasing, which is useful when plotting the first few coordinates. When \code{\link{makeFV}} is applied to a test set, \code{Xf} are coordinates of the feature vectors of the test set in the same space as those of the training set, and then \code{kmat} = \code{Xf} \%*\% \code{Xf of training data}.
  }
\item{transfmat}{square matrix for transforming kmat to \code{Xf}.}
}

\references{
Raymaekers J., Rousseeuw P.J., Hubert M. (2021). Class maps for visualizing classification results. \emph{Technometrics}, appeared online. \doi{10.1080/00401706.2021.1927849}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J., Hubert, M.
}

\seealso{
\code{\link{makeKernel}}
}

\examples{
library(e1071)
set.seed(1); X <- matrix(rnorm(200 * 2), ncol = 2)
X[1:100, ] <- X[1:100, ] + 2
X[101:150, ] <- X[101:150, ] - 2
y <- as.factor(c(rep("blue", 150), rep("red", 50)))
cols <- c("deepskyblue3", "red")
plot(X, col = cols[as.numeric(y)], pch = 19)
# We now fit an SVM with radial basis kernel to the data:
svmfit <- svm(y~.,  data = data.frame(X = X, y = y),  scale = FALSE,
             kernel = "radial", cost = 10, gamma = 1, probability = TRUE)
Kxx <- makeKernel(X, svfit = svmfit)
outFV <- makeFV(Kxx)
Xf <- outFV$Xf # The data matrix in this feature space.
dim(Xf) # The feature vectors are high dimensional.
# The inner products of Xf match the kernel matrix:
max(abs(as.vector(Kxx - crossprod(t(Xf), t(Xf))))) # 3.005374e-13 # tiny, OK
range(rowSums(Xf^2)) # all points in Xf lie on the unit sphere.
pairs(Xf[, 1:5], col = cols[as.numeric(y)])
# In some of these we see spherical effects, e.g.
plot(Xf[, 1], Xf[, 5], col = cols[as.numeric(y)], pch = 19)
# The data look more separable here than in the original
# two-dimensional space.

# For more examples, we refer to the vignette:
\dontrun{
vignette("Support_vector_machine_examples")
}
}

