\name{calcTemperatedLogDen.default}
\alias{calcTemperatedLogDen.default}
\title{calcTemperatedLogDen default}
\description{Rescale Log-Densities by given Temperatures}
\usage{
\method{calcTemperatedLogDen}{default}(x,
temp,
...,
refDen = apply(x,
2,
quantile,
probs = 0.9))
}
\arguments{
  \item{x}{numeric matrix (nStep x nResComp): logDensities for each component at each step}
  \item{temp}{numeric vector (nResComp): Temperature, i.e. cost reduction factor}
  \item{\dots}{
}
  \item{refDen}{reference density: differences to this value will be scaled and readded}
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
\value{numeric matrix (nStep x nResComp), rescaled logDen}

\author{Thomas Wutzler}



\seealso{\code{\link{twDEMC}}

\code{\link{sumLogDenComp}}}

