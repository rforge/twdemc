\name{initZtwDEMCSub.twDEMC}
\alias{initZtwDEMCSub.twDEMC}
\title{initZtwDEMCSub twDEMC}
\description{generates an appropriate initial sample of parameter vectors for twDEMC from subsampling a previous result}
\usage{
\method{initZtwDEMCSub}{twDEMC}(vtwdemc,
...,
vars = colnames(vtwdemc$parms),
nChainPop = getNChainsPop(vtwdemc),
nPop = getNPops(vtwdemc))
}
\arguments{
  \item{vtwdemc}{the twDEMC list to subsample}
  \item{\dots}{
}
  \item{vars}{which variables to keep}
  \item{nChainPop}{
}
  \item{nPop}{
}
}



\author{Thomas Wutzler}



\seealso{\code{\link{initZtwDEMCNormal}}}

