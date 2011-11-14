\name{calcDEMCTempDiffLogDen2Init}
\alias{calcDEMCTempDiffLogDen2Init}
\title{calcDEMCTempDiffLogDen2Init}
\description{Estimate scalar temperatures to obtain acceptance rate }
\usage{
calcDEMCTempDiffLogDen2Init(resLogDen, ..., doConstrainNeg = TRUE)
}
\arguments{
  \item{resLogDen}{result of \code{\link{twCalcLogDenPar}}: list with components  logDen and logDenComp}
  \item{\dots}{further arguments to \code{\link{calcDEMCTempDiffLogDenConst}}}
  \item{doConstrainNeg}{different default value}
}

\value{named numeric vector of estimated Temperatures per data stream. }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





