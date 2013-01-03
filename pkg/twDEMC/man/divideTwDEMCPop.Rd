\name{divideTwDEMCPop}
\alias{divideTwDEMCPop}
\title{divideTwDEMCPop}
\description{Splits a specific populations of x into subpopulations for given subSpaces. }
\usage{divideTwDEMCPop(pop, subSpaces, isCheckBorderConsistency = TRUE)}
\arguments{
  \item{pop}{numeric array (nStep, nParm, nChain): an object of class aTwDEMCPops returned by \code{\link{twDEMCBlockInt}}}
  \item{subSpaces}{a list (nPop): each entry listing subspaces for the given population with entries lowerParBounds, upperParBounds, and splits}
  \item{isCheckBorderConsistency}{if FALSE the check for border consistency is omitted.}
}



\author{Thomas Wutzler}




\examples{
data(den2dCorEx)
x <- den2dCorEx$mcBulk
newPops <- divideTwDEMCPop(x$pops[[1]], den2dCorEx$subspaces0[[1]]$spaces )
length(newPops)
str(newPops[[1]])
#lapply(newPops,"[[","splits")
}
