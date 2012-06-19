\name{getNGen.twDEMCPops}
\alias{getNGen.twDEMCPops}
\title{getNGen twDEMCPops}
\description{Extract the number of completed generations in res}
\usage{
\method{getNGen}{twDEMCPops}(x,
...)
}
\arguments{
  \item{x}{object of class twDEMCPops}
  \item{\dots}{
}
}
\details{the number of generations corresponds to the steps after time 0 (sample 1).
Hence the sample of size 2 and thinning 1 describes one generation (one step forward).
A sample of size 2 of thinning 5 encompasses times 0 and 5, i.e. 5 generations.
see }


\author{Thomas Wutzler}



\seealso{\code{\link{iSample2time}}
\code{\link{time2iSample}}
\code{\link{getNPops.twDEMCPops}}
\code{\link{getNSamples.twDEMCPops}}
\code{\link{getNSamplesSpace.twDEMCPops}}
\code{\link{getNChains.twDEMCPops}}
\code{\link{getNChainsPop.twDEMCPops}}
\code{\link{getNParms.twDEMCPops}}
\code{\link{getNSpaces.twDEMCPops}}
\code{\link{getSpacesPop.twDEMCPops}}
\code{\link{subset.twDEMCPops}}
,\code{\link{twDEMCBlockInt}}}
\examples{
data(twdemcEx1)
getNGen(twdemcEx1)
getNSamples(twdemcEx1)
getNSamplesPop(twdemcEx1)
getNSamplesSpace(twdemcEx1)
twdemcEx1$thin
getNPops(twdemcEx1)
getNChains(twdemcEx1)
getNChainsPop(twdemcEx1)
getNParms(twdemcEx1)
getNBlocks(twdemcEx1)
}
