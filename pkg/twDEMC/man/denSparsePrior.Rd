\name{denSparsePrior}
\alias{denSparsePrior}
\title{denSparsePrior}
\description{Same as denSpace but returning two components: prior and misfit}
\usage{denSparsePrior(theta, twTwoDenEx = twTwoDenEx1, theta0 = twTwoDenEx$thetaTrue, 
    thetaPrior = NULL, invCovarTheta = NULL, intermediate = list(), 
    ...)}
\arguments{
  \item{theta}{
}
  \item{twTwoDenEx}{
}
  \item{theta0}{
}
  \item{thetaPrior}{the prior estimate of the parameters}
  \item{invCovarTheta}{the inverse of the Covariance of the prior parameter estimates}
  \item{intermediate}{specify a mutable environement to store intermediate results that are the same across different densities for given theta}
  \item{\dots}{
}
}



\author{Thomas Wutzler}





