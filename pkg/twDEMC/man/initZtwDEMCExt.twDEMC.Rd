\name{initZtwDEMCExt.twDEMC}
\alias{initZtwDEMCExt.twDEMC}
\title{initZtwDEMCExt twDEMC}
\description{subsampling and extending twDEMC with new variables}
\usage{\method{initZtwDEMCExt}{twDEMC}(vtwdemc, ..., thetaPrior, covarTheta, nChainPop = getNChainsPop(vtwdemc), 
    nPop = getNPops(vtwdemc))}
\arguments{
  \item{vtwdemc}{the twDEMC list to subsample}
  \item{\dots}{e.g. m0}
  \item{thetaPrior}{numeric vector: mu of multivariate gaussian prior distribtuion}
  \item{covarTheta}{numeric vector: sigma of multivariate gaussian prior distrbituion}
  \item{nChainPop}{
}
  \item{nPop}{
}
}



\author{Thomas Wutzler}



\seealso{\code{\link{initZtwDEMCNormal}}}

