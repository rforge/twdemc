\name{subsetTail.twDEMCPops}
\alias{subsetTail.twDEMCPops}
\title{subsetTail twDEMCPops}
\description{discards the first part of all the chains}
\usage{
\method{subsetTail}{twDEMCPops}(x,
pKeep = 0.5,
...)
}
\arguments{
  \item{x}{object of class twDEMCPops}
  \item{pKeep}{the percentage of the samples to keep}
  \item{\dots}{
}
}



\author{Thomas Wutzler}




\examples{
data(twdemcEx1)
res <- subsetTail(twdemcEx1)
plot( as.mcmc.list(res), smooth=FALSE )
getNSamples(twdemcEx1)
getNSamples(res)

# special case: no rows in one pop of twdemcEx1
twdemcEx1$pops[[2]] <- .subsetTwDEMCPop
}
