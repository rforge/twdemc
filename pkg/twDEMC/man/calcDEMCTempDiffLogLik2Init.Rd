\name{calcDEMCTempDiffLogLik2Init}
\alias{calcDEMCTempDiffLogLik2Init}
\title{calcDEMCTempDiffLogLik2Init}
\description{Estimate scalar temperatures to obtain acceptance rate }
\usage{
calcDEMCTempDiffLogLik2Init(resLogLik, ..., doConstrainNeg = TRUE)
}
\arguments{
  \item{resLogLik}{result of \code{\link{twCalcLogLikPar}}: list with components  logLik and resFLogLik}
  \item{\dots}{further arguments to \code{\link{calcDEMCTempDiffLogLikConst}}}
  \item{doConstrainNeg}{different default value}
}

\value{named numeric vector of estimated Temperatures per data stream. }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





