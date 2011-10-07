\name{calcDEMCTempDiffLogDenInit}
\alias{calcDEMCTempDiffLogDenInit}
\title{calcDEMCTempDiffLogDenInit}
\description{Estimate Temperatures for different data streams to obtain acceptance rate }
\usage{
calcDEMCTempDiffLogDenInit(resLogDen, ...)
}
\arguments{
  \item{resLogDen}{result of \code{\link{twCalcLogDenPar}}: list with components  logDen and resFLogDen}
  \item{\dots}{further arguments to \code{\link{calcDEMCTempDiffLogDen}}}
}

\value{named numeric vector of estimated Temperatures per data stream. }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





