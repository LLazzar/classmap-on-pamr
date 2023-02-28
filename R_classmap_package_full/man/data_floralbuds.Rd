\name{data_floralbuds}
\alias{data_floralbuds}
\docType{data}
\title{
Floral buds data
}
\description{
This data on floral pear bud detection was first described by Wouters et al.
The goal is to classify the instances into buds, branches, scales and support.
The numeric vectors resulted from a multispectral vision sensor and describe the scanned images.
}
\usage{data("data_floralbuds")}

\format{
  A data frame with 550 observations on the following 7 variables.
  \describe{
    \item{\code{X1}}{numeric vector}
    \item{\code{X2}}{numeric vector}
    \item{\code{X3}}{numeric vector}
    \item{\code{X4}}{numeric vector}
    \item{\code{X5}}{numeric vector}
    \item{\code{X6}}{numeric vector}
    \item{\code{y}}{a factor with levels \code{branch} \code{bud} \code{scales} \code{support}}
  }
}


\source{
Wouters, N., De Ketelaere, B., Deckers, T. De Baerdemaeker, J., Saeys, W. (2015).
Multispectral detection of floral buds for automated thinning of pear.
\emph{Comput. Electron. Agric.} 113, C, 93–103. <doi:10.1016/j.compag.2015.01.015>
}

\examples{
data("data_floralbuds")
str(data_floralbuds)
summary(data_floralbuds)

# The data are used in:
\dontrun{
vignette("Discriminant_analysis_examples")
vignette("Neural_net_examples")
}
}
