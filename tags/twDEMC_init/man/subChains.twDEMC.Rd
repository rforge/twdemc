\name{subChains.twDEMC}
\alias{subChains.twDEMC}
\title{subChains twDEMC}
\description{Condenses an twDEMC List to the chains iChains e.g \code{1:4}.}
\usage{
subChains.twDEMC(x, iChains = NULL, iPops = NULL, nPop = NULL, ..., doKeepBatchCall = FALSE)
}
\arguments{
  \item{x}{
}
  \item{iChains}{
}
  \item{iPops}{
}
  \item{nPop}{number of populations in x, if not specified then taken from ncol(x$temp)}
  \item{\dots}{
}
  \item{doKeepBatchCall}{wheter to retain the batch call attribute of x}
}
\details{There are several methods to transform or subset the results of an \code{\link{twDEMCInt}} run. \itemize{
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
\item{ plotting a subset of the chains and cases: \code{\link{plotThinned.mcmc.list}}  } 
\item{ transforming parameters \code{\link{transOrigPopt.mcmc.list}}  } 
}

Alternatively to specification of iChains, one can specify a vector of populations and the total number of populations.}
\value{a list of class twDEMC (see \code{\link{twDEMCInt}})}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twDEMCInt}}}

