\name{replaceNonFiniteDiffLogDens}
\alias{replaceNonFiniteDiffLogDens}
\title{replaceNonFiniteDiffLogDens}
\description{For each component replace NAs by sample of others and non-Finite values by minimum of others  }
\usage{
replaceNonFiniteDiffLogDens(diffLogDen, doConstrainNeg = FALSE)
}
\arguments{
  \item{diffLogDen}{array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}}
  \item{doConstrainNeg}{if given, density of accepted jumps (positive) is constrained to 0}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





