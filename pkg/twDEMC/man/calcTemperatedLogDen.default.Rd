\name{calcTemperatedLogDen.default}
\alias{calcTemperatedLogDen.default}
\title{calcTemperatedLogDen default}
\description{Rescale Log-Densities by given Temperatures}
\usage{\method{calcTemperatedLogDen}{default}(x, temp, ...)}\arguments{
  \item{x}{numeric matrix (nStep x nResComp): logDensities for each component at each step}
  \item{temp}{numeric vector (nResComp): Temperature, i.e. cost reduction factor}
  \item{\dots}{
}
}
\details{There are several function to help with temperating log-Densities
\itemize{
\item{ Rescale Log-Densities by given Temperatures: this function  } 
\item{ Rescale Log-Densities within chains \code{\link{calcTemperatedLogDenChains.array}} }
\item{ Generating exponential temperature series: \code{\link{calcDEMCTemp}}  }
\item{ Calculating a base temperature from stream temperatures: \code{\link{calcBaseTemp}}  }
\item{ Calculating a base temperature from cost values: \code{\link{calcBaseTempSk}}  }
\item{ Calculating a Stream temperatures from base temperature: \code{\link{calcStreamTemp}}  }
}}
\value{numeric matrix (nStep x nResComp), rescaled logDen}

\author{Thomas Wutzler}



\seealso{\code{\link{twDEMC}}}

