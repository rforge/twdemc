\name{stackChainsPop.twDEMCPops}
\alias{stackChainsPop.twDEMCPops}
\title{stackChainsPop twDEMCPops}
\description{Combine MarkovChains of each population of a twDEMC to a single chain.}
\usage{\method{stackChainsPop}{twDEMCPops}(x, mergeMethod = "stack", ...)}\arguments{
  \item{x}{
}
  \item{mergeMethod}{
}
  \item{\dots}{
}
}



\author{Thomas Wutzler}




\examples{
data(twdemcEx1)
getNChainsPop(twdemcEx1)	# four chains within population
getNSamples(twdemcEx1)		# with 26 samples

res <- stackChainsPop(twdemcEx1)
getNChainsPop(res)			# only 1 chains
getNSamples(res)			# but with 26*4=104 samples
str(concatPops(res))
}
