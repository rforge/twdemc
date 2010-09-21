\name{thin.twDEMC}
\alias{thin.twDEMC}
\title{thin twDEMC}
\description{Reduces the rows of an twDEMC object (list returned by \code{\link{twDEMCInt}}) to correspond to a thinning of \code{newThin}.}
\usage{
thin.twDEMC(x, newThin = x$thin, start = 1, end = NULL, ..., doKeepBatchCall = FALSE)
}
\arguments{
  \item{x}{the twDEMC list to thin}
  \item{newThin}{the target thinning factor, must be positive multiple of vMcpl$thin}
  \item{start}{the start time of the chain}
  \item{end}{the maximum end time of the chains}
  \item{\dots}{
}
  \item{doKeepBatchCall}{wheter to retain the batch call attribute of x}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{subChains.twDEMC}}}

