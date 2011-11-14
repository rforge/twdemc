\name{calcDEMCTempDiffLogDen2}
\alias{calcDEMCTempDiffLogDen2}
\title{calcDEMCTempDiffLogDen2}
\description{Estimate scalar temperature vector obtain given acceptance rate and scaling temperatures to the same magnitude }
\usage{
calcDEMCTempDiffLogDen2(diffLogDen, pTarget = 0.18, TFix = numeric(0), Tmax = .Machine$double.xmax, rFracMin = 1/4, doConstrainNeg = FALSE)
}
\arguments{
  \item{diffLogDen}{array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}}
  \item{pTarget}{minimum acceptance rate of component}
  \item{TFix}{named numeric vector: components whose Temperate is to be fixed}
  \item{Tmax}{maximum temperature to return}
  \item{rFracMin}{fraction of Density to which data-streams with low diff-logDen are scaled to by temperatures to}
  \item{doConstrainNeg}{if given, density of accepted jumps (positive) is constrained to 0}
}

\value{numeric vector of result components 
yielding acceptance rates closest to pTarget}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





