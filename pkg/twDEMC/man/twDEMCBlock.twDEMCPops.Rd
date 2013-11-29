\name{twDEMCBlock.twDEMCPops}
\alias{twDEMCBlock.twDEMCPops}
\title{twDEMCBlock twDEMCPops}
\description{initialize \code{\link{twDEMCBlockInt}} by former run and append results to former run}
\usage{\method{twDEMCBlock}{twDEMCPops}(x, ..., TEnd = numeric(0), TEnd0, TStart = numeric(0), 
    TStart0, TFix, nObs, doRecordProposals = FALSE, extendRun = TRUE)}
\arguments{
  \item{x}{list of class twDEMCPops, result of \code{\link{twDEMCBlockInt}}}
  \item{\dots}{further arguments to \code{\link{twDEMCBlockInt}}}
  \item{TEnd}{numeric vector (nResComp) of end temperatures for each result component, if not given stays at current temperature, if in addition TSpec is specified, it has no effect}
  \item{TEnd0}{numeric scalar: alternative way of specifying end tempertaure: by baseTemperature for all data streams}
  \item{TStart}{numeric vector (nResComp): specifing starting temperature, if neigher TSpec, TStart, or TStart0 are given stays at the current temperature}
  \item{TStart0}{numeric scalar: specifing starting temperature by baseTemperature for all data streams}
  \item{TFix}{see \code{\link{twDEMCSACont}}, argument \code{ctrlT$TFix}. Must be given if specifying TEnd0 or TStart0}
  \item{nObs}{see \code{\link{twDEMCSACont}}, argument must be given if specifying TEnd0 or TStart0}
  \item{doRecordProposals}{if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y}
  \item{extendRun}{if set to FALSE then only the new samples are returned}
}
\details{pops, dInfos, and blocks are reused from x or overwritten by arguments

\describe{\item{TSpec}{ 
If TSpec is not explicitly given as an array for \code{\link{twDEMCBlockInt}} it is constructed as follows:
, If TEnd is given use it to explicitely define end Temperature for all data streams.
, If TEnd0 is given calculate TEnd for all data streams by scaling for number of observations (see \code{\link{calcStreamTemp}}).
, If TStart0 is given calcualte TCurr for all data streams 
, else use current Temperature of streams in x as starting temperature
}}}


\author{Thomas Wutzler}




\examples{
data(twdemcEx1) 		# previous run of twDEMCBlock
class(twdemcEx1)
twdemcEx1$thin			# thinning interval
(nGen0 <- getNGen(twdemcEx1))		# number of generations

# extend by 16 generations
nGen <- 16
#mtrace(twDEMCBlock.twDEMCPops)
res <- twDEMCBlock( twdemcEx1, nGen=nGen )

identical( nGen0+nGen, getNGen(res) )
plot( as.mcmc.list(res), smooth=FALSE )
}
