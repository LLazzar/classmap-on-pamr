\name{data_bookReviews}
\alias{data_bookReviews}
\docType{data}
\title{
Amazon book reviews data
}
\description{
This is a subset of the data used in the paper, which
was assembled by Prettenhofer and Stein (2010). It contains 1000 reviews
of books on Amazon, of which 500 were selected from the original training data
and 500 from the test data.

The full dataset has been used for a variety of things, including
classification using svm. The subset was chosen small enough to keep the computation
time low, while still containing the examples in the paper.
}

\usage{data("data_bookReviews")}

\format{
  A data frame with 1000 observations on the following 2 variables.
  \describe{
    \item{\code{review}}{the review in text format (character)}
    \item{\code{sentiment}}{factor indicating the sentiment of the review: negative (1) or positive (2)}
  }
}


\source{
Prettenhofer, P., Stein, B. (2010). Cross-language text classification using
structural correspondence learning. \emph{Proceedings of the 48th annual
meeting of the association for computational linguistics},  1118-1127.
}

\examples{
data(data_bookReviews)
# Example review:
data_bookReviews[5, 1]

# The data are used in:
\dontrun{
vignette("Support_vector_machine_examples")
}
}
