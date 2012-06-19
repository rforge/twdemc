\name{subset.twDEMCPops}
\alias{subset.twDEMCPops}
\title{subset twDEMCPops}
\description{Condenses an twDEMCPops result object to the cases boKeep.}
\usage{
\method{subset}{twDEMCPops}(x,
boKeep,
...,
iPops = seq_along(x$pops),
dropShortPops = FALSE)
}
\arguments{
  \item{x}{twDEMCPops object}
  \item{boKeep}{either logical vector or numeric vector of indices of cases to keep}
  \item{\dots}{
}
  \item{iPops}{integer vector: only these populations are subset, others are kept}
  \item{dropShortPops}{if set to TRUE, pops in iPops with less samples than to what \code{boKeep} refers to are dropped}
}
\details{There are several methods access properties a result of an \code{\link{twDEMCBlockInt}}, i.e.
an object of class \code{twDEMCPops} \itemize{
\item{ number of generations: \code{\link{getNGen.twDEMCPops}}  } 
\item{ number of samples (only one sample each thinning inteval): \code{\link{getNSamples.twDEMCPops}}  } 
\item{ number of chains: \code{\link{getNChains.twDEMCPops}}  } 
\item{ number of populations: \code{\link{getNPops.twDEMCPops}}  } 
\item{ number of chains per population: \code{\link{getNChainsPop.twDEMCPops}}  } 
\item{ number of parameters: \code{\link{getNParms.twDEMCPops}}  } 
\item{ thinning interval: \code{res$thin}  } 
\item{ space replicate that the poplation belongs to: \code{\link{spacesPop.twDEMCPops}}  } 
\item{ number of space replicates: \code{\link{getNSpaces.twDEMCPops}}  } 
\item{ number of blocks: \code{\link{getNBlocks.twDEMCPops}}  } 
\item{ parameter bounds: \code{\link{getParBoundsPop.twDEMCPops}}  } 
\item{ parameter bounds: \code{\link{getCurrentTemp.twDEMCPops}}  } 
}

There are several methods to transform or subset the results of an \code{\link{twDEMCBlockInt}} run. \itemize{
\item{ transforming to array representation across populations of type \code{twDEMC}: code{\link{concatPops.twDEMCPops}}  }
\item{ select sub-populations: \code{\link{subPops.twDEMCPops}}  } 
\item{ select subset of cases: \code{\link{subset.twDEMCPops}} (this function) }
\item{ make populations the same length: \code{\link{squeeze.twDEMCPops}}  }
\item{ stack all the results of all chains to one big matrix: \code{\link{stackChains.twDEMCPops}}  } 
\item{ thin all the chains: \code{\link{thin.twDEMCPops}}  } 
\item{ combine several twDEMCPops results to a bigger set of populations: \code{\link{combineTwDEMCPops}}  }
\item{ combine populations for subspaces to bigger populations: \code{\link{stackPopsInSpace.twDEMCPops}}  }
\item{ combine MarkovChains of each population of a twDEMC to a single chain: \code{\link{stackChainsPop.twDEMCPops}}  }
}

There are several methods utilize the functions of the coda package. \itemize{
\item{ convert an twDEMCPops to a coda mcmc.list \code{\link{as.mcmc.list.twDEMCPops}}  } 
\item{ applying a function to all of the chains: \code{\link{mcmcListApply}}  }
\item{ stack all the results of all chains to one big matrix: \code{\link{stackChains.mcmc.list}}  } 
\item{ transforming parameters \code{\link{transOrigPopt.mcmc.list}}  } 
}

components \code{thin,Y,nGenBurnin} are kept, but may be meaningless after subsetting.}


\author{Thomas Wutzler}



\seealso{\code{\link{twDEMCBlockInt}}}
\examples{
if( FALSE){
	tmp <- subset(res,1:3)
	str(tmp)
	tmp <- subset(res,1:10)  # should produce an error
	#mtrace(subset.twDEMCPops)
	tmp <- subset(res,1:10,dropShortPops=TRUE)  
} 
}
