\name{calcDEMCTempDiffLogDenConst}
\alias{calcDEMCTempDiffLogDenConst}
\title{calcDEMCTempDiffLogDenConst}
\description{Estimate scalar temperature to obtain given acceptance rate }
\usage{
calcDEMCTempDiffLogDenConst(diffLogDen, pTarget = 0.18, TFix = numeric(0), Tmax = .Machine$double.xmax)
}
\arguments{
  \item{diffLogDen}{array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}}
  \item{pTarget}{minimum acceptance rate of component}
  \item{TFix}{named numeric vector: components whose Temperate is to be fixed}
  \item{Tmax}{maximum temperature to return}
}

\value{numeric scalar Temperature between in [1,Tmax] 
yielding acceptance rates closest to pTarget}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





