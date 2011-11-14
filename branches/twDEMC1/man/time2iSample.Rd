\name{time2iSample}
\alias{time2iSample}
\title{time2iSample}
\description{Convert time to sample number}
\usage{
time2iSample(time, thin = 1, match = c(round = "round", floor = "floor", ceiling = "ceiling", none = "none"))
}
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

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{iSample2time}}
\code{\link{getNGen.twDEMC}}
,\code{\link{subChains.twDEMC}}
,\code{\link{twDEMCInt}}}

