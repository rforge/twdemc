\name{calcTemperatedLogDen.twDEMCPops}
\alias{calcTemperatedLogDen.twDEMCPops}
\title{calcTemperatedLogDen twDEMCPops}
\description{Rescale Log-Densities by given Temperatures}
\usage{\method{calcTemperatedLogDen}{twDEMCPops}(x, temp = getCurrentTemp(x), 
    ...)}\arguments{
  \item{x}{object of class twDEMCPops}
  \item{temp}{numeric vector (nResComp): Temperature, i.e. cost reduction factor}
  \item{\dots}{
}
}

\value{array (nCases * nComp) of stacked rescaled logDensities}

\author{Thomas Wutzler}



\seealso{\code{\link{calcTemperatedLogDen.default}}

\code{\link{calcTemperatedLogDenChains.twDEMC}}}

