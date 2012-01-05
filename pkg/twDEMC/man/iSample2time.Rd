\name{iSample2time}
\alias{iSample2time}
\title{iSample2time}
\description{Convert sample number to time given the thinning interval}
\usage{
iSample2time(iSample,
thin = 1)
}
\arguments{
  \item{iSample}{integer vector: indices of the samples in recorded states}
  \item{thin}{integer scalar: thinning interval}
}
\details{Sample 1 corresponds to time zero
After first generation, the time sample 2 corresponds to thin}


\author{Thomas Wutzler}



\seealso{\code{\link{time2iSample}}
\code{\link{getNGen.twDEMC}}
,\code{\link{subChains.twDEMC}}
,\code{\link{twDEMCBlockInt}}}
\examples{
iSample <- 1:13
structure( iSample2time(iSample) , names=iSample )
structure( iSample2time(iSample, thin=10) , names=iSample )
structure( iSample2time(iSample, thin=4) , names=iSample )
structure( iSample2time(iSample, thin=8) , names=iSample )
}
