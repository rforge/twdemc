\name{twRunDivideTwDEMCBatch}
\alias{twRunDivideTwDEMCBatch}
\title{twRunDivideTwDEMCBatch}
\description{wrapper to divideTwDEMC compliant to runCluster.R}
\usage{
twRunDivideTwDEMCBatch(argsDivideTwDEMCBatch, ..., prevResRunCluster = NULL, restartFilename = NULL)
}
\arguments{
  \item{argsDivideTwDEMCBatch}{Arguments passed to DivideTwDEMCBatch -> DivideTwDEMC -> twDEMCInt -> fLogDen.
It is updated by \dots.
After update it must contain entries Zinit and fLogDen
It is further searched for entries nPops, and argsFLogDen, and logDenCompX. The latter are initialized to defaults  \code{1,list(),character(0)} respectively if not found.   }
  \item{\dots}{further arguments passed to DivideTwDEMCBatch -> DivideTwDEMC -> twDEMCInt -> fLogDen}
  \item{prevResRunCluster}{results of call to twRunDEMC, argument required to be called from runCluster.R}
  \item{restartFilename}{name of the file to store restart information, argument required to be called from runCluster.R}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





