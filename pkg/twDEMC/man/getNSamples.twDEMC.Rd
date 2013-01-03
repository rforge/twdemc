\name{getNSamples.twDEMC}
\alias{getNSamples.twDEMC}
\title{getNSamples twDEMC}
\description{Extracts the number of s}
\usage{\method{getNSamples}{twDEMC}(res, ...)}
\arguments{
  \item{res}{object of class twDEMC}
  \item{\dots}{
}
}
\details{There is only one sample per thinning interval of length \code{res$thin}.}


\author{Thomas Wutzler}



\seealso{\code{\link{getNGen.twDEMC}}
\code{\link{subChains.twDEMC}}
,\code{\link{twDEMCBlockInt}}}
\examples{FALSE}
