\name{combinePops.twDEMC}
\alias{combinePops.twDEMC}
\title{combinePops twDEMC}
\description{Combine several populations to one big population.}
\usage{
combinePops.twDEMC(x, ..., doKeepBatchCall = FALSE)
}
\arguments{
  \item{x}{first twDEMC object}
  \item{\dots}{more twDEMC objects}
  \item{doKeepBatchCall}{wheter to retain the batch call attribute of x}
}
\details{All arguments must be of type twDEMC and of same thinning interval and length.}


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{subChains.twDEMC}}}

