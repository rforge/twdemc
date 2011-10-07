\name{subset.twDEMC}
\alias{subset.twDEMC}
\title{subset twDEMC}
\description{Condenses an twDEMC result object to the cases boKeep.}
\usage{
\method{subset}{twDEMC}(x, boKeep, ...)
}
\arguments{
  \item{x}{twDEMC object}
  \item{boKeep}{either logical vector or numeric vector of indices of cases to keep}
  \item{\dots}{
}
}
\details{If last row is dropped then all components of resFLogLikX are set to -Inf
This ensures that it LogLik for last state is recalculated if object is used again in twDEMC

components \code{thin,Y,nGenBurnin} are kept, but may be meaningless after subsetting.}
\value{list of class twDEMC with subset of cases in parsm, rLogLik, pAccept, and temp}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{subChains.twDEMC}}}

