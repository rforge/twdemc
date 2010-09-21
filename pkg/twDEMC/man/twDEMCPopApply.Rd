\name{twDEMCPopApply}
\alias{twDEMCPopApply}
\title{twDEMCPopApply}
\description{Applying a function across all chains of one population for each case.}
\usage{
twDEMCPopApply(x, nPops, FUN, ...)
}
\arguments{
  \item{x}{a matrix with columns chains or array with last dimension chain}
  \item{nPops}{number of populations}
  \item{FUN}{function to apply to population submatrix}
  \item{\dots}{further arguemtns to FUN}
}

\value{array with last dimenstion correponding to population}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





