\name{calcTemperatedLogDen.default}
\alias{calcTemperatedLogDen.default}
\title{calcTemperatedLogDen default}
\description{Rescale Log-Densities by given Temperatures}
\usage{\method{calcTemperatedLogDen}{default}(x, temp, ..., maxLogDen = apply(x, 
    2, max))}
\arguments{
  \item{x}{numeric matrix (nStep x nResComp): logDensities for each component at each step}
  \item{temp}{numeric vector (nResComp): Temperature, i.e. cost reduction factor}
  \item{\dots}{
}
  \item{maxLogDen}{best achievable log-Density}
}
\details{There are several function to help with temperating log-Densities
\itemize{
\item{ Rescale Log-Densities by given Temperatures: this function  } 
\item{ Generating exponential T-series: \code{\link{calcDEMCTemp}}  }
}
Others functions that deal with TProp are deprecated. 


Only the difference of the log-Density to the best achievable Log-Density 
is important for comparison.
Hence, first the best logDen is subtracted, then the difference is multiplied by given temperatue.}


\author{Thomas Wutzler}



\seealso{\code{\link{twDEMC}}}

