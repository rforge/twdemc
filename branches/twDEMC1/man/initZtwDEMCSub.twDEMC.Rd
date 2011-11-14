\name{initZtwDEMCSub.twDEMC}
\alias{initZtwDEMCSub.twDEMC}
\title{initZtwDEMCSub twDEMC}
\description{generates an appropriate initial sample of parameter vectors for twDEMC from subsampling a previous result}
\usage{
\method{initZtwDEMCSub}{twDEMC}(vtwdemc, vars = rownames(vtwdemc$parms), nChains = ncol(vtwdemc$rLogDen), nPops = ncol(vtwdemc$temp), ...)
}
\arguments{
  \item{vtwdemc}{the twDEMC list to subsample}
  \item{vars}{which variables to keep}
  \item{nChains}{
}
  \item{nPops}{
}
  \item{\dots}{
}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}}

