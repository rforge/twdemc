\name{stackChains.twDEMC}
\alias{stackChains.twDEMC}
\title{stackChains twDEMC}
\description{Combine MarkovChains of a twDEMC to a matrix. }
\usage{\method{stackChains}{twDEMC}(x, omitBurnin = FALSE, ...)}\arguments{
  \item{x}{
}
  \item{omitBurnin}{if TRUE, then burnin of each chain is omitted before stacking}
  \item{\dots}{
}
}

\value{Matrix with first attributes$nBlock columns the logDensity logDen and the remaining columns the variables.}

\author{Thomas Wutzler}



\seealso{\code{\link{subChains.twDEMC}}}

