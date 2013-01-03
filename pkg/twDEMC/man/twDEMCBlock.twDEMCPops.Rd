\name{twDEMCBlock.twDEMCPops}
\alias{twDEMCBlock.twDEMCPops}
\title{twDEMCBlock twDEMCPops}
\description{initialize \code{\link{twDEMCBlockInt}} by former run and append results to former run}
\usage{\method{twDEMCBlock}{twDEMCPops}(x, ..., TEnd = numeric(0), doRecordProposals = FALSE, 
    extendRun = TRUE)}
\arguments{
  \item{x}{list of class twDEMC, result of \code{\link{twDEMCBlockInt}}}
  \item{\dots}{further arguments to \code{\link{twDEMCBlockInt}}}
  \item{TEnd}{numeric vector (nResComp) of end temperatures for each result component, if not given stays at current temperature, if in addition TSpec is specified, it has no effect}
  \item{doRecordProposals}{if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y}
  \item{extendRun}{if set to FALSE then only the new samples are returned}
}
\details{pops, dInfos, and blocks are reused from x or overwritten by arguments}


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
