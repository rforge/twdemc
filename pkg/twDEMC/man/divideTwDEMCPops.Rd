\name{divideTwDEMCPops}
\alias{divideTwDEMCPops}
\title{divideTwDEMCPops}
\description{splits all populations of x into subpopulations for given subSpaces }
\usage{divideTwDEMCPops(x, subSpacesPop)}
\arguments{
  \item{x}{an object of class aTwDEMCPops returned by \code{\link{twDEMCBlockInt}}}
  \item{subSpacesPop}{a list (nPop): each entry listing subspaces for the given population}
}



\author{Thomas Wutzler}




\examples{
data(den2dCorEx)
#x <- den2dCorEx$mcSubspaces0
x <- den2dCorEx$mcBulk
#sapply( x$pops, "[[", "spaceInd")
#subSpaces <- do.call( c, lapply(den2dCorEx$subspaces0,"[[","spaces"))
subSpacesPop <- den2dCorEx$subspaces0
#mtrace(divideTwDEMCPops)
xNew <- divideTwDEMCPops(x, subSpacesPop )
sapply(den2dCorEx$subspaces0, function(subSpacesPop){ length(subSpacesPop$spaces) } )
getNPops(xNew)
getSpacesPop(xNew)
str(xNew$pops[[9]])
#lapply(xNew$pops,"[[","splits")
}
