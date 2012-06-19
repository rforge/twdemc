\name{as.mcmc.list.twDEMCPops}
\alias{as.mcmc.list.twDEMCPops}
\title{as mcmc list twDEMCPops}
\description{Converts list of type twDEMCPops (result of \code{\link{twDEMCBlockInt}}) to coda's \code{mcmc.list}. }
\usage{
\method{as.mcmc}{list.twDEMCPops}(x,
...,
useThinning = TRUE,
minPopLength = NULL)
}
\arguments{
  \item{x}{the output of \code{\link{twDEMCBlockInt}}) run}
  \item{\dots}{
}
  \item{useThinning}{if TRUE thinning is used to make populations the same length, if FALSE they are cut to shortest population}
  \item{minPopLength}{integer: if given, shorter populations are dropped}
}



\author{Thomas Wutzler}





