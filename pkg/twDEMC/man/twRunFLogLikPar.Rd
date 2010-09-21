\name{twRunFLogLikPar}
\alias{twRunFLogLikPar}
\title{twRunFLogLikPar}
\description{Wrapper for \code{\link{twCalcLogLikPar}} to accept arguments provided by runCluster.R}
\usage{
twRunFLogLikPar(..., prevResRunCluster = NULL, restartFilename = NULL)
}
\arguments{
  \item{\dots}{further arguments passed \code{\link{twCalcLogLikPar}}}
  \item{prevResRunCluster}{results of call to twRunDEMC, argument required to be called from runCluster.R}
  \item{restartFilename}{name of the file to store restart information, argument required to be called from runCluster.R}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





