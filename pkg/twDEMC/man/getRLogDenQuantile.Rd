\name{getRLogDenQuantile}
\alias{getRLogDenQuantile}
\title{getRLogDenQuantile}
\description{Quantile of logDensity below which models are significantly different from the best model, i.e. parameterization}
\usage{
getRLogDenQuantile(stackedSample, maxLogDen = max(stackedSample[, 1]), df = ncol(stackedSample) - 1, perc = 0.95)
}
\arguments{
  \item{stackedSample}{numeric matrix: first column logDensity, other columns free parameters, see \code{\link{stackChains.twDEMC}}}
  \item{maxLogDen}{maximum logDen Density}
  \item{df}{degress of freedom: number of fitted parameters}
  \item{perc}{percentile of significance}
}
\details{See Hilborn97 for explanation of Density ratio test for nested models.}
\value{numeric scalar: minimum LogDensity below which models are significantly different }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{checkConvergenceGelman}}
\code{\link{twDEMCInt}}}

