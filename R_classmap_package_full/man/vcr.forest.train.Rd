\name{vcr.forest.train}
\alias{vcr.forest.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare for visualization of a random forest classification on training data
}
\description{
Produces output for the purpose of constructing graphical displays such as the \code{\link{classmap}} and \code{\link{silplot}}. The user first needs to train a random forest on the data by \code{\link[randomForest:randomForest]{randomForest::randomForest}}.
This then serves as an argument to  \code{\link{vcr.forest.train}}.
}

\usage{
vcr.forest.train(X, y, trainfit, type = list(),
                 k = 5, stand = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A rectangular matrix or data frame, where the
             columns (variables) may be of mixed type.
}
\item{y}{factor with the given class labels.
             It is crucial that \code{X} and \code{y} are exactly the same
             as in the call to \code{\link[randomForest:randomForest]{randomForest::randomForest}}.
             \code{y} is allowed to contain \code{NA}'s.}
\item{trainfit}{the output of a \code{\link[randomForest:randomForest]{randomForest::randomForest}} training run.}
\item{k}{the number of nearest neighbors used in the
             farness computation.}
\item{type}{list for specifying some (or all) of the types of the
    variables (columns) in \code{X}, used for computing the dissimilarity matrix, as in \code{\link[cluster:daisy]{cluster::daisy}}.  The list may contain the following components: \code{"ordratio"} (ratio scaled variables to be treated as ordinal variables), \code{"logratio"} (ratio scaled variables that
    must be logarithmically transformed), \code{"asymm"} (asymmetric
    binary) and \code{"symm"} (symmetric binary variables).  Each
    component's value is a vector, containing the names or the numbers
    of the corresponding columns of \code{X}.
    Variables not mentioned in the \code{type} list are interpreted as
    usual (see argument \code{X}).}
\item{stand}{whether or not to standardize numerical (interval scaled) variables by their range as in the original \code{\link[cluster:daisy]{cluster::daisy}} code for the farness computation. Defaults to \code{TRUE}.}
}


\value{
A list with components: \cr
  \item{X}{The data used to train the forest.}
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
  \item{trainfit}{The trained random forest which was given as an input to this function.}
}
\references{
Raymaekers J., Rousseeuw P.J.(2021). Silhouettes and quasi residual plots for neural nets and tree-based classifiers. \href{https://arxiv.org/abs/2106.08814}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.forest.newdata}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}
}

\examples{
library(randomForest)
data("data_instagram")
traindata <- data_instagram[which(data_instagram$dataType == "train"), -13]
set.seed(71) # randomForest is not deterministic
rfout <- randomForest(y~., data = traindata, keep.forest = TRUE)
mytype <- list(symm = c(1, 5, 7, 8)) # These 4 columns are
# (symmetric) binary variables. The variables that are not
# listed are interval-scaled by default.
x_train <- traindata[, -12]
y_train <- traindata[, 12]
# Prepare for visualization:
vcrtrain <- vcr.forest.train(X = x_train, y = y_train,
                            trainfit = rfout, type = mytype)
confmat.vcr(vcrtrain)
stackedplot(vcrtrain, classCols = c(4, 2))
silplot(vcrtrain, classCols = c(4, 2))
classmap(vcrtrain, "genuine", classCols = c(4, 2))
classmap(vcrtrain, "fake", classCols = c(4, 2))

# For more examples, we refer to the vignette:
\dontrun{
vignette("Random_forest_examples")
}
}
