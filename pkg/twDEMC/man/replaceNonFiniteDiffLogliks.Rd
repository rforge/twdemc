\name{replaceNonFiniteDiffLogliks}
\alias{replaceNonFiniteDiffLogliks}
\title{replaceNonFiniteDiffLogliks}
\description{For each component replace NAs by sample of others and non-Finite values by minimum of others  }
\usage{
replaceNonFiniteDiffLogliks(diffLogLik, doConstrainNeg = FALSE)
}
\arguments{
  \item{diffLogLik}{array( streams x steps) Lp-La see \code{\link{getDiffLogLik.twDEMCProps}}}
  \item{doConstrainNeg}{if given, likelihood of accepted jumps (positive) is constrained to 0}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





