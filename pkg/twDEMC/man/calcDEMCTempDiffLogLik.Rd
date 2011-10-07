\name{calcDEMCTempDiffLogLik}
\alias{calcDEMCTempDiffLogLik}
\title{calcDEMCTempDiffLogLik}
\description{Estimate temperatures for different data streams to obtain given acceptance rate }
\usage{
calcDEMCTempDiffLogLik(diffLogLik, pTarget = 0.18, TFix = numeric(0))
}
\arguments{
  \item{diffLogLik}{array( streams x steps x  chains) Lp-La see \code{\link{getDiffLogLik.twDEMCProps}}}
  \item{pTarget}{overall acceptance rate}
  \item{TFix}{named numeric vector: components whose Temperate is to be fixed}
}

\value{numeric vector of Temperatures (fLogLik_Component x chains)}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





