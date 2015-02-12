\name{calcStreamTemp}
\alias{calcStreamTemp}
\title{calcStreamTemp}
\description{scale base temperature to given number of observations.}
\usage{calcStreamTemp(tempBase, nObs, TFix = rep(NA_real_, length(nObs)), 
    iFixTemp = which(is.finite(TFix)))}\arguments{
  \item{tempBase}{numeric scalar: The temperature, i.e. variance inflation factor, for a single observation.}
  \item{nObs}{integer vector (nResultComp): number of observations}
  \item{TFix}{numeric vector (nResultComp): fixed temperature for components, non-finite for those with varying temperature}
  \item{iFixTemp}{integer vector: index of result components, which Temperature is fixed}
}

\value{numeric vector: The Temperatures for the observation streams.}

\author{Thomas Wutzler}



\seealso{\link{calcBaseTemp}, \link{twDEMC}

\code{\link{calcTemperatedLogDen.default}}}

