\name{calcDEMCTempDiffLogLik3Init}
\alias{calcDEMCTempDiffLogLik3Init}
\title{calcDEMCTempDiffLogLik3Init}
\description{Estimate scalar temperatures to obtain acceptance rate }
\usage{
calcDEMCTempDiffLogLik3Init(resLogLik, ..., doConstrainNeg = TRUE)
}
\arguments{
  \item{resLogLik}{result of \code{\link{twCalcLogLikPar}}: list with components  logLik and resFLogLik}
  \item{\dots}{further arguments to \code{\link{calcDEMCTempDiffLogLikConst}}}
  \item{doConstrainNeg}{different default value}
}

\value{named numeric vector of estimated Temperatures per data stream. }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





