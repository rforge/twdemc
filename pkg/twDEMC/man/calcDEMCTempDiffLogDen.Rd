\name{calcDEMCTempDiffLogDen}
\alias{calcDEMCTempDiffLogDen}
\title{calcDEMCTempDiffLogDen}
\description{Estimate temperatures for different data streams to obtain given acceptance rate }
\usage{
calcDEMCTempDiffLogDen(diffLogDen, pTarget = 0.18, TFix = numeric(0))
}
\arguments{
  \item{diffLogDen}{array( streams x steps x  chains) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}}
  \item{pTarget}{overall acceptance rate}
  \item{TFix}{named numeric vector: components whose Temperate is to be fixed}
}

\value{numeric vector of Temperatures (fLogDen_Component x chains)}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





