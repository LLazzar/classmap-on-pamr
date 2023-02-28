\name{vcr.rpart.train}
\alias{vcr.rpart.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare for visualization of an rpart classification on training data.
}
\description{
Produces output for the purpose of constructing graphical displays such as the \code{\link{classmap}}. The user first needs to train a
classification tree on the data by \code{\link[rpart:rpart]{rpart::rpart}}.
This then serves as an argument to  \code{\link{vcr.rpart.train}}.
}

\usage{
vcr.rpart.train(X, y, trainfit, type = list(),
                k = 5, stand = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A rectangular matrix or data frame, where the
             columns (variables) may be of mixed type and
             may contain \code{NA}'s.
}
\item{y}{factor with the given class labels.
             It is crucial that \code{X} and \code{y} are exactly the same
             as in the call to \code{\link[rpart:rpart]{rpart::rpart}}.
             \code{y} is allowed to contain \code{NA}'s.}
\item{k}{the number of nearest neighbors used in the
             farness computation.}
\item{trainfit}{the output of an \code{\link[rpart:rpart]{rpart::rpart}} training cycle.}
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
  \item{X}{The input data \code{X}. Keep??}
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
 \item{trainfit}{the trainfit used to build the VCR object.}
}
\references{
Raymaekers J., Rousseeuw P.J.(2021). Silhouettes and quasi residual plots for neural nets and tree-based classifiers. \href{https://arxiv.org/abs/2106.08814}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.rpart.newdata}}, \code{\link{classmap}}, \code{\link{silplot}}, \code{\link{stackedplot}}
}

\examples{
library(rpart)
data("data_titanic")
traindata <- data_titanic[which(data_titanic$dataType == "train"), -13]
str(traindata); table(traindata$y)
set.seed(123) # rpart is not deterministic
rpart.out <- rpart(y ~ Pclass + Sex + SibSp +
                    Parch + Fare + Embarked,
                  data = traindata, method = 'class', model = TRUE)
y_train <- traindata[, 12]
x_train <- traindata[, -12]
mytype <- list(nominal = c("Name", "Sex", "Ticket", "Cabin", "Embarked"), ordratio = c("Pclass"))
# These are 5 nominal columns, and one ordinal.
# The variables not listed are by default interval-scaled.
vcrtrain <- vcr.rpart.train(x_train, y_train, rpart.out, mytype)
confmat.vcr(vcrtrain)
silplot(vcrtrain, classCols = c(2, 4))
classmap(vcrtrain, "casualty", classCols = c(2, 4))
classmap(vcrtrain, "survived", classCols = c(2, 4))

# For more examples, we refer to the vignette:
\dontrun{
vignette("Rpart_examples")
}
}
