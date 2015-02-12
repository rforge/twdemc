\name{calcBaseTempSk}
\alias{calcBaseTempSk}
\title{calcBaseTempSk}
\description{Calculate the mean variance factor from given cost and observation numbers.        }
\usage{calcBaseTempSk(Sk, nObs = 1, TFix = rep(NA_real_, ncol(Sk)), 
    iFixTemp = which(is.finite(TFix)), iNonFixTemp = which(!is.finite(TFix)), 
    isVerbose = FALSE)}\arguments{
  \item{Sk}{numeric matrix (nCases, nResultComp): misfit (=-2*logDensity) per result components}
  \item{nObs}{number of observations per result component}
  \item{TFix}{numeric vector (nResultComp): fixed temperature for components, non-finite for those with varying temperature}
  \item{iFixTemp}{integer vector: index of result components, which Temperature is fixed}
  \item{iNonFixTemp}{integer vector: index of result components, which Temperature is fixed}
  \item{isVerbose}{set to TRUE to report mean base temperatures per stream}
}

\value{numeric scalar: The temperature, i.e. variance inflation factor, for a single observation.}

\author{Thomas Wutzler}



\seealso{\code{\link{calcTemperatedLogDen.default}}}

