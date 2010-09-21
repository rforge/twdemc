\name{initZtwDEMCExt.twDEMC}
\alias{initZtwDEMCExt.twDEMC}
\title{initZtwDEMCExt twDEMC}
\description{subsampling and extending twDEMC with new variables}
\usage{
initZtwDEMCExt.twDEMC(vtwdemc, thetaPrior, covarTheta, nChains = ncol(vtwdemc$rLogLik), nPops = ncol(vtwdemc$temp), ...)
}
\arguments{
  \item{vtwdemc}{the twDEMC list to subsample}
  \item{thetaPrior}{numeric vector: mu of multivariate gaussian prior distribtuion}
  \item{covarTheta}{numeric vector: sigma of multivariate gaussian prior distrbituion}
  \item{nChains}{
}
  \item{nPops}{
}
  \item{\dots}{e.g. m0}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}}

