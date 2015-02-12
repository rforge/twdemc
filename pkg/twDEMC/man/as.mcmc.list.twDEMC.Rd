\name{as.mcmc.list.twDEMC}
\alias{as.mcmc.list.twDEMC}
\title{as mcmc list twDEMC}
\description{Converts list of type twDEMC (result of \code{\link{twDEMC}}) to coda's \code{mcmc.list}. }
\usage{\method{as.mcmc}{list.twDEMC}(x, maxLength = NULL, thin = 1, start = 1, 
    ...)}\arguments{
  \item{x}{the output of \code{\link{twDEMC}}) run}
  \item{maxLength}{maximum length of the chains, thin will be increased accordingly}
  \item{thin}{thinning interval}
  \item{start}{starting generation (for neglecting the burnin)}
  \item{\dots}{
}
}

\value{a \code{\link[coda]{mcmc.list}}}

\author{Thomas Wutzler}



\seealso{\code{\link{subChains.twDEMC}}}

