\name{initZtwDEMCExt.matrix}
\alias{initZtwDEMCExt.matrix}
\title{initZtwDEMCExt matrix}
\description{subsampling and extending twDEMC with new variables}
\usage{
\method{initZtwDEMCExt}{matrix}(Zinit1, thetaPrior, covarTheta, nChains = 4, nPops = 1, m0 = calcM0twDEMC(length(thetaPrior), nPops, nChains), ...)
}
\arguments{
  \item{Zinit1}{the twDEMC list to subsample}
  \item{thetaPrior}{numeric vector: mu of multivariate gaussian prior distribtuion}
  \item{covarTheta}{numeric vector: sigma of multivariate gaussian prior distrbituion}
  \item{nChains}{
}
  \item{nPops}{
}
  \item{m0}{number of required cases }
  \item{\dots}{
}
}
\details{the variables in thetaPrior that are part of vtwdemc are subsampled
the other variables are drawn from prior distribution
assuming no correlations between variables present and absent in vtwdemc and }


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}}

