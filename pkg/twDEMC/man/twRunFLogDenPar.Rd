\name{twRunFLogDenPar}
\alias{twRunFLogDenPar}
\title{twRunFLogDenPar}
\description{Wrapper for \code{\link{twCalcLogDenPar}} to accept arguments provided by runCluster.R}
\usage{
twRunFLogDenPar(..., prevResRunCluster = NULL, restartFilename = NULL)
}
\arguments{
  \item{\dots}{further arguments passed \code{\link{twCalcLogDenPar}}}
  \item{prevResRunCluster}{results of call to twRunDEMC, argument required to be called from runCluster.R}
  \item{restartFilename}{name of the file to store restart information, argument required to be called from runCluster.R}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





