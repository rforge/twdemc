\name{calcDEMCTempDiffLogDenConstInit}
\alias{calcDEMCTempDiffLogDenConstInit}
\title{calcDEMCTempDiffLogDenConstInit}
\description{Estimate scalar temperatures to obtain acceptance rate }
\usage{
calcDEMCTempDiffLogDenConstInit(resLogDen, ...)
}
\arguments{
  \item{resLogDen}{result of \code{\link{twCalcLogDenPar}}: list with components  logDen and logDenComp}
  \item{\dots}{further arguments to \code{\link{calcDEMCTempDiffLogDenConst}}}
}

\value{named numeric vector of estimated Temperatures per data stream. }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





