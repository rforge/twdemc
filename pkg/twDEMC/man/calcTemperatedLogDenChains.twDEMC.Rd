\name{calcTemperatedLogDenChains.twDEMC}
\alias{calcTemperatedLogDenChains.twDEMC}
\title{calcTemperatedLogDenChains twDEMC}
\description{Rescale Log-Densities by given Temperatures within chains}
\usage{
\method{calcTemperatedLogDenChains}{twDEMC}(x,
temp = getCurrentTemp(x),
...)
}
\arguments{
  \item{x}{object with entry resLogDen, a numeric matrix (nStep x nResComp x nChain), e.g. twDEMC or pop in twDEMCPops}
  \item{temp}{numeric vector (nResComp): Temperature, i.e. cost reduction factor}
  \item{\dots}{
}
}



\author{Thomas Wutzler}





