\name{stackChains.twDEMC}
\alias{stackChains.twDEMC}
\title{stackChains twDEMC}
\description{Combine MarkovChains of a twDEMC to a matrix. }
\usage{
\method{stackChains}{twDEMC}(x, omitBurnin = FALSE, ...)
}
\arguments{
  \item{x}{
}
  \item{omitBurnin}{if TRUE, then burnin of each chain is omitted before stacking}
  \item{\dots}{
}
}

\value{Matrix with first column the logDensity rLogDen and the remaining columns the variables.}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{subChains.twDEMC}}}

