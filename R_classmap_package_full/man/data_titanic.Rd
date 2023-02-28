\name{data_titanic}
\alias{data_titanic}
\docType{data}
\title{
Titanic data
}
\description{

This dataset contains information on 1309 passengers of the RMS Titanic.
The goal is to predict survival based on 11 characteristics such as the travel class, age and sex of the passengers.

The original data source is
https://www.kaggle.com/c/titanic/data

The data is split up in a training data consisting of 891 observations and a test data of 418 observations. The response in the test set was obtained by combining information from other data files, and has been verified by submitting it as a `prediction' to kaggle and getting perfect marks.
}

\usage{data("data_titanic")}

\format{
  A data frame with 1309 observations on the following variables.
  \describe{
  \item{PassengerId}{a unique identified for each passenger.}
\item{Pclass}{travel class of the passenger.}
\item{Name}{name of the passenger.}
\item{Sex}{sex of the passenger.}
\item{Age}{age of the passenger.}
\item{SibSp}{number of siblings and spouses traveling with the passenger.}
\item{Parch}{number of parents and children traveling with the passenger.}
\item{Ticket}{Ticket number of the passenger.}
\item{Fare}{fare paid for the ticket.}
\item{Cabin}{cabin number of the passenger.}
\item{Embarked}{Port of embarkation. Takes the values C (Cherbourg), Q (Queenstown) and S (Southampton).}
\item{y}{factor indicating casualty or survivor.}
  \item{dataType}{vector taking the values ``train'' or ``test'' indicating whether the observation belongs to the training or the test data.}
  }
}


\source{
https://www.kaggle.com/c/titanic/data
}

\examples{
data("data_titanic")
traindata <- data_titanic[which(data_titanic$dataType == "train"), -13]
testdata <- data_titanic[which(data_titanic$dataType == "test"), -13]
str(traindata)
table(traindata$y)

# The data are used in:
\dontrun{
vignette("Rpart_examples")
}
}
