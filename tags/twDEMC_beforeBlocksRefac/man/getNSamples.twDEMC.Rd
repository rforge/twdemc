\name{getNSamples.twDEMC}
\alias{getNSamples.twDEMC}
\title{getNSamples twDEMC}
\description{Extracts the number of s}
\usage{
\method{getNSamples}{twDEMC}(res, ...)
}
\arguments{
  \item{res}{object of class twDEMC}
  \item{\dots}{
}
}
\details{There is only one sample per thinning interval of length \code{res$thin}.}
\value{integer, number of samples in twDEMC}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{getNGen.twDEMC}}
\code{\link{subChains.twDEMC}}
,\code{\link{twDEMCInt}}}

