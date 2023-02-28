\name{data_instagram}
\alias{data_instagram}
\docType{data}
\title{
Instagram data
}
\description{

This dataset contains information on fake (spam) accounts on Instagram.
The original source is https://www.kaggle.com/free4ever1/instagram-fake-spammer-genuine-accounts by Bardiya Bakhshandeh.

The data contains information on 696 Instagram accounts.
For each account, 11 variables were recorded describing its characteristics.
The goal is to detect fake instagram accounts, which are used for spamming.
}

\usage{data("data_instagram")}

\format{
  A data frame with 696 observations on the following variables.
  \describe{
  \item{profile.pic}{binary, indicates whether profile has picture.}
\item{nums.length.username}{ratio of number of numerical chars in username to its length.}
\item{fullname.words}{number of words in full name.}
\item{nums.length.fullname}{ratio of number of numerical characters in full name to its length.}
\item{name..username}{binary, indicates whether the name and username of the profile
are the same.}
\item{description.length}{length of the description/biography of the profile (in number of characters).}
\item{external.URL}{binary, indicates whether profile has external url.}
\item{private}{binary, indicates whether profile is private or not.}
\item{X.posts}{number of posts made by profile.}
\item{X.followers}{number of followers.}
\item{X.follows}{numbers of follows.}
\item{y}{whether profile is fake or not.}
  \item{dataType}{vector taking the values ``train'' or ``test'' indicating
  whether the observation belongs to the training or the test data.}
  }
}


\source{
https://www.kaggle.com/free4ever1/instagram-fake-spammer-genuine-accounts
}

\examples{
data(data_instagram)
str(data_instagram)

# The data are used in:
\dontrun{
vignette("Random_forest_examples")
}
}
