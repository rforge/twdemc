\name{denBoth}
\alias{denBoth}
\title{denBoth}
\description{Example of using two different logDensity functions: comparison of combining denSparce and denRich into one function.}
\usage{denBoth(theta, twTwoDenEx = twTwoDenEx1, weights = c(1, 1, 1), 
    theta0 = twTwoDenEx$thetaTrue, thetaPrior = NULL, invCovarTheta = NULL, 
    ...)}
\arguments{
  \item{theta}{
}
  \item{twTwoDenEx}{
}
  \item{weights}{
}
  \item{theta0}{
}
  \item{thetaPrior}{the prior estimate of the parameters}
  \item{invCovarTheta}{the inverse of the Covariance of the prior parameter estimates}
  \item{\dots}{
}
}



\author{Thomas Wutzler}




\examples{
data(twTwoDenEx1)
denBoth(twTwoDenEx1$thetaTrue)
denBoth(twTwoDenEx1$thetaTrue, thresholdCovar = 0.3)
}
