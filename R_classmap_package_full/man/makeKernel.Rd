\name{makeKernel}
\alias{makeKernel}
\title{
Compute kernel matrix
}
\description{
Computes kernel value or kernel matrix, where the kernel type is extracted from an svm trained by \code{\link[e1071:svm]{e1071::svm}}.
}
\usage{
makeKernel(X1, X2 = NULL, svfit)
}

\arguments{

\item{X1}{first matrix (or vector) of coordinates.}
\item{X2}{if not \code{NULL}, second data matrix or vector. If NULL, \code{X2} is assumed equal to \code{X1}.}
\item{svfit}{output from \code{\link[e1071:svm]{e1071::svm}}}.
}


\value{the kernel matrix, of dimensions \code{nrow(X1)} by \code{nrow(X2)}. When both \code{X1} and \code{X2} are vectors, the result is a single number.
}

\references{
Raymaekers J., Rousseeuw P.J., Hubert M. (2021). Class maps for visualizing classification results. \emph{Technometrics}, appeared online. \doi{10.1080/00401706.2021.1927849}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{makeFV}}
}

\examples{

library(e1071)
set.seed(1); X <- matrix(rnorm(200 * 2), ncol = 2)
X[1:100, ] <- X[1:100, ] + 2
X[101:150, ] <- X[101:150, ] - 2
y <- as.factor(c(rep("blue", 150), rep("red", 50))) # two classes
# We now fit an SVM with radial basis kernel to the data:
set.seed(1) # to make the result of svm() reproducible.
svmfit <- svm(y~.,  data = data.frame(X = X, y = y),  scale = FALSE,
             kernel = "radial", cost = 10, gamma = 1, probability = TRUE)
Kxx <- makeKernel(X, svfit = svmfit)
# The result is a square kernel matrix:
dim(Kxx) # 200 200
Kxx[1:5, 1:5]

# For more examples, we refer to the vignette:
\dontrun{
vignette("Support_vector_machine_examples")
}
}

