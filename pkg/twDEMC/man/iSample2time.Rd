\name{iSample2time}
\alias{iSample2time}
\title{iSample2time}
\description{Convert sample number to time given the thinning interval}
\usage{
iSample2time(iSample, thin = 1)
}
\arguments{
  \item{iSample}{integer vector: indices of the samples in recorded states}
  \item{thin}{integer scalar: thinning interval}
}
\details{Sample 1 corresponds to time zero
After first generation, the time sample 2 corresponds to thin}


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{time2iSample}}
\code{\link{getNGen.twDEMC}}
,\code{\link{subChains.twDEMC}}
,\code{\link{twDEMCInt}}}
\examples{
time = 0:20
structure( time2iSample(time), names=time)
structure( iSample2time(1:6, thin=4) , names=1:6 )	# to show the times corresponding to sample
structure( time2iSample(time,thin=4,match="floor"), names=time)
structure( time2iSample(time,thin=4,match="ceiling"), names=time)
structure( time2iSample(time,thin=4), names=time)	# round may vary for times exactly between two samples
structure( time2iSample(time,thin=4,match="none"), names=time)	# here we do not get indices
}
