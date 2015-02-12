\name{denRichPrior}
\alias{denRichPrior}
\title{denRichPrior}
\description{Example of using two different logDensity functions: density of data-rich observations}
\usage{denRichPrior(theta, twTwoDenEx, theta0 = twTwoDenEx$thetaTrue, 
    thetaPrior = NULL, invCovarTheta = NULL, intermediate = list(), 
    ...)}\arguments{
  \item{theta}{
}
  \item{twTwoDenEx}{
}
  \item{theta0}{
}
  \item{thetaPrior}{the prior estimate of the parameters}
  \item{invCovarTheta}{the inverse of the Covariance of the prior parameter estimates, or alternatively the diag of Covariance matrix (sigma_i)}
  \item{intermediate}{specify a mutable environement to store intermediate results that are the same across different densities for given theta}
  \item{\dots}{
}
}



\author{Thomas Wutzler}



\seealso{\code{\link{denSparsePrior}}}

