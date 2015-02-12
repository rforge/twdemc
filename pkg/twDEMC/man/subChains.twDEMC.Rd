\name{subChains.twDEMC}
\alias{subChains.twDEMC}
\title{subChains twDEMC}
\description{Condenses an twDEMC List to the chains iChains e.g \code{1:4}.}
\usage{\method{subChains}{twDEMC}(x, iChains = NULL, iPops = NULL, nPop = getNPops(x), 
    ..., doKeepBatchCall = FALSE)}\arguments{
  \item{x}{the twDEMC object to select chains from}
  \item{iChains}{integer vector: the chains to select}
  \item{iPops}{integer vector: the populations, whose chains should be selected}
  \item{nPop}{number of populations in x}
  \item{\dots}{not used}
  \item{doKeepBatchCall}{wheter to retain the batch call attribute of x}
}
\details{Alternatively to specification of \code{iChains}, one can specify a vector of populations and the total number of populations.

\describe{ \item{The twDEMC class}{
The \code{twDEMC} class represents MCMC results with all chains having the same length. 
This is usually a result by a call to \code{\link{concatPops.twDEMCPops}}.
It holds the same results as the \code{twDEMCPops} class, except that the entries in \code{pops} list are concatenated by chain.  
\itemize{
\item parms: array (nStep, nParm, nChain) of parameter samples
\item temp: matrix (nStep, nResComp ) of temperature, one collumn for each logDensity components
\item parms: array (nStep, nBlock, nChain) acceptance rate
\item resLogDen: array (nStep, nResComp, nChain) results of logDensities
\item logDen: array (nStep, nBlock, nChain) sum across logDensity components for each block
\item Y: array (lastUnthinnedSteps, nParm+nBlock+nResComp, nChain) of parameters, acceptance per block and result components of latest unthinned steps 
}

There are several methods access properties, i.e.
an object of class \code{twDEMC} \itemize{
\item{ number of generations: \code{\link{getNGen.twDEMC}}  } 
\item{ number of samples (only one sample each thinning inteval): \code{\link{getNSamples.twDEMC}}  } 
\item{ number of chains: \code{\link{getNChains.twDEMC}}  } 
\item{ number of populations: \code{\link{getNPops.twDEMC}}  } 
\item{ number of chains per population: \code{\link{getNChainsPop.twDEMC}}  } 
\item{ number of parameters: \code{\link{getNParms.twDEMC}}  } 
\item{ thinning interval: \code{res$thin}  } 
}

There are several methods to transform or subset: \itemize{
\item{ select chains or sub-populations: this method  } 
\item{ thin all the chains: \code{\link{thin.twDEMC}}  } 
\item{ select subset of cases: \code{\link{subset.twDEMC}}  }
\item{ combine several twDEMC results to a bigger set of populations \code{\link{combinePops.twDEMC}}  }
\item{ stack all the results of all chains to one big matrix \code{\link{stackChains.twDEMC}}  } 
}

There are several methods utilize the functions of the coda package. \itemize{
\item{ convert an twDEMC to a coda mcmc.list \code{\link{as.mcmc.list.twDEMC}}  } 
\item{ applying a function to all of the chains: \code{\link{mcmcListApply}}  }
\item{ stack all the results of all chains to one big matrix: \code{\link{stackChains.mcmc.list}}  } 
\item{ transforming parameters \code{\link{transOrigPopt.mcmc.list}}  } 
}
}}}
\value{A list of class twDEMC (see \code{\link{twDEMCBlockInt}}).}

\author{Thomas Wutzler}



\seealso{\code{\link{twDEMC}}}

