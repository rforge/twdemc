\name{stackPopsInSpace.twDEMCPops}
\alias{stackPopsInSpace.twDEMCPops}
\title{stackPopsInSpace twDEMCPops}
\description{Combine populations for subspaces to bigger populations }
\usage{\method{stackPopsInSpace}{twDEMCPops}(x, ..., spacePop = getSpacesPop(x))}
\arguments{
  \item{x}{
}
  \item{\dots}{arguments passed \code{\link{combineTwDEMCPops}} such as \code{mergeMethod=stack/slice/random}}
  \item{spacePop}{
}
}



\author{Thomas Wutzler}



\seealso{\code{\link{stackPops.twDEMCPops}}
\code{\link{combineTwDEMCPops}}
\code{\link{stackChains.twDEMCPops}}
\code{\link{subset.twDEMCPops}}}
\examples{
data(den2dCorEx)
getNSamples(den2dCorEx$mcBulk)
res <- stackPopsInSpace( den2dCorEx$mcSubspaces0 )
getNSamples(res)	# lost a few samples in sorting chains to subspaces
#mtrace(concatPops.twDEMCPops)
plot( as.mcmc.list(den2dCorEx$mcBulk), smooth=FALSE ) # original before splitting into subspaces
plot( as.mcmc.list(res), smooth=FALSE )		# stacked populations of subspaces
plot( as.mcmc.list(stackPopsInSpace( den2dCorEx$mcSubspaces0,mergeMethod="stack" )), smooth=FALSE )
plot( as.mcmc.list(stackPopsInSpace( den2dCorEx$mcSubspaces0,mergeMethod="slice" )), smooth=FALSE )
}
