\name{concatPops.twDEMCPops}
\alias{concatPops.twDEMCPops}
\title{concatPops twDEMCPops}
\description{Concatenates all the chains of all subpopulations to array across all chains as class \code{twDEMC}.}
\usage{\method{concatPops}{twDEMCPops}(x,
 ..., isUsingThinning = TRUE, minPopLength = NULL)}\arguments{
  \item{x}{the twDEMCPops object to transform}
  \item{\dots}{not used}
  \item{isUsingThinning}{if TRUE (defaul) thinning is used to make populations the same length (the minimum across populations' length), if FALSE they are cut to shortest population}
  \item{minPopLength}{integer scalar: if specified, populations with less samples than length.out are dropped from the results}
}
\details{In the twDEMCPops object \code{x}, the information on results is scattered in a list of populations 
(result component \code{pop} described in \code{link{twDEMCBlockInt}}).
This function makes all chains the same length, and combines the populations by appending all the chains in a big array.
All other entry besides \code{pops} is retained from the original twDEMCPops object \code{x}.}
\value{An object of class \code{twDEMC} (see \code{\link{subChains.twDEMC}})}

\author{Thomas Wutzler}



\seealso{\code{\link{subset.twDEMCPops}}
,\code{\link{subChains.twDEMC}}
,\code{\link{twDEMCBlockInt}}}

