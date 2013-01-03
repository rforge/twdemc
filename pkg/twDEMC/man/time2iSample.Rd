\name{time2iSample}
\alias{time2iSample}
\title{time2iSample}
\description{Convert time to sample number}
\usage{time2iSample(time, thin = 1, match = c(round = "round", floor = "floor", 
    ceiling = "ceiling", none = "none"))}
\arguments{
  \item{time}{numeric: indices of the samples in recorded states}
  \item{thin}{integer scalar: thinning interval}
  \item{match}{mode of matching times between samples\describe{
\item{round}{closest sample}
\item{floor}{sample before time}
\item{ceiling}{sample after time}
\item{none}{returns a fractional, use if you sure that time corresponds to actual sample indices}
}}
}

\value{index of the neares sample}

\author{Thomas Wutzler}



\seealso{\code{\link{iSample2time}}
\code{\link{getNGen.twDEMC}}
,\code{\link{subChains.twDEMC}}
,\code{\link{twDEMCBlockInt}}}
\examples{
time = 0:20
structure( time2iSample(time), names=time)
structure( iSample2time(1:6, thin=4) , names=1:6 )	# to show the times corresponding to sample
structure( time2iSample(time,thin=4,match="floor"), names=time)
structure( time2iSample(time,thin=4,match="ceiling"), names=time)
structure( time2iSample(time,thin=4), names=time)	# round may vary for times exactly between two samples
structure( time2iSample(time,thin=4,match="none"), names=time)	# here we do not get indices
}
