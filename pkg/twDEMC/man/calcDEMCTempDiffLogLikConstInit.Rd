\name{calcDEMCTempDiffLogLikConstInit}
\alias{calcDEMCTempDiffLogLikConstInit}
\title{calcDEMCTempDiffLogLikConstInit}
\description{Estimate scalar temperatures to obtain acceptance rate }
\usage{
calcDEMCTempDiffLogLikConstInit(resLogLik, ...)
}
\arguments{
  \item{resLogLik}{result of \code{\link{twCalcLogLikPar}}: list with components  logLik and resFLogLik}
  \item{\dots}{further arguments to \code{\link{calcDEMCTempDiffLogLikConst}}}
}

\value{named numeric vector of estimated Temperatures per data stream. }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





