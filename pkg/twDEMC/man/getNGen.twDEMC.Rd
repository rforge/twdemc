\name{getNGen.twDEMC}
\alias{getNGen.twDEMC}
\title{getNGen twDEMC}
\description{Extract the number of completed generations in res}
\usage{\method{getNGen}{twDEMC}(res, ...)}
\arguments{
  \item{res}{object of class twDEMC}
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
\code{\link{getNPops.twDEMC}}
\code{\link{getNSamples.twDEMC}}
\code{\link{getNChains.twDEMC}}
\code{\link{getNChainsPop.twDEMC}}
\code{\link{getNParms.twDEMC}}
\code{\link{subChains.twDEMC}}
,\code{\link{twDEMCBlockInt}}}
\examples{
data(twdemcEx1)
getNGen(twdemcEx1)
getNSamples(twdemcEx1)
twdemcEx1$thin
getNPops(twdemcEx1)
getNChains(twdemcEx1)
getNChainsPop(twdemcEx1)
getNParms(twdemcEx1)
getNBlocks(twdemcEx1)
}
