\name{calcDEMCTempDiffLogLik3}
\alias{calcDEMCTempDiffLogLik3}
\title{calcDEMCTempDiffLogLik3}
\description{Estimate scalar temperature vector obtain given acceptance rate and scaling temperatures to the same magnitude }
\usage{
calcDEMCTempDiffLogLik3(diffLogLik, pTarget = 0.18, TFix = numeric(0), Tmax = .Machine$double.xmax, rFracMin = 1/4, doConstrainNeg = FALSE)
}
\arguments{
  \item{diffLogLik}{array( streams x steps) Lp-La see \code{\link{getDiffLogLik.twDEMCProps}}}
  \item{pTarget}{minimum acceptance rate of component}
  \item{TFix}{named numeric vector: components whose Temperate is to be fixed}
  \item{Tmax}{maximum temperature to return}
  \item{rFracMin}{fraction of Likelihood to which data-streams with low diff-logLik are scaled to by temperatures to}
  \item{doConstrainNeg}{if given, likelihood of accepted jumps (positive) is constrained to 0}
}

\value{numeric vector of result components Temperatures
yielding acceptance rates closest to pTarget
attribute pAcceptTVar giving the calculated acceptance rate for Temperature dependent step}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





