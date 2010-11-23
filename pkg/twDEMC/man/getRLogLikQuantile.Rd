\name{getRLogLikQuantile}
\alias{getRLogLikQuantile}
\title{getRLogLikQuantile}
\description{Quantile of logLikelihood below which models are significantly different from the best model, i.e. parameterization}
\usage{
getRLogLikQuantile(stackedSample, maxLogLik = max(stackedSample[, 1]), df = ncol(stackedSample) - 1, perc = 0.95)
}
\arguments{
  \item{stackedSample}{numeric matrix: first column log-Likelihood, other columns free parameters, see \code{\link{stackChains.twDEMC}}}
  \item{maxLogLik}{maximum logLik Likelihood}
  \item{df}{degress of freedom: number of fitted parameters}
  \item{perc}{percentile of significance}
}
\details{See Hilborn97 for explanation of Likelihood ratio test for nested models.}
\value{numeric scalar: minimum Log-Likelihood below which models are significantly different }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{checkConvergenceGelman}}
\code{\link{twDEMCInt}}}

