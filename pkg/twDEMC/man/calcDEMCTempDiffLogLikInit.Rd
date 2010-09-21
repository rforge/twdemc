\name{calcDEMCTempDiffLogLikInit}
\alias{calcDEMCTempDiffLogLikInit}
\title{calcDEMCTempDiffLogLikInit}
\description{Estimate Temperatures for different data streams to obtain acceptance rate }
\usage{
calcDEMCTempDiffLogLikInit(resLogLik, ...)
}
\arguments{
  \item{resLogLik}{result of \code{\link{twCalcLogLikPar}}: list with components  logLik and resFLogLik}
  \item{\dots}{further arguments to \code{\link{calcDEMCTempDiffLogLik}}}
}

\value{named numeric vector of estimated Temperatures per data stream. }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





