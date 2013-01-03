\name{combinePops.twDEMC}
\alias{combinePops.twDEMC}
\title{combinePops twDEMC}
\description{Combine several populations to one big population consiting of more chains.}
\usage{\method{combinePops}{twDEMC}(x, ..., doKeepBatchCall = FALSE)}
\arguments{
  \item{x}{first twDEMC object}
  \item{\dots}{more twDEMC objects}
  \item{doKeepBatchCall}{wheter to retain the batch call attribute of x}
}
\details{All arguments must be of type twDEMC and of same thinning interval and length.
Temperature is taken from the first population - it assumes that all populations have the same Temperature in their steps}


\author{Thomas Wutzler}



\seealso{\code{\link{subChains.twDEMC}}}
\examples{FALSE}
