\name{twRunDEMC}
\alias{twRunDEMC}
\title{twRunDEMC}
\description{Run a twDEMC with first replacing nonfinite Likelihoods in Zinits last row.}
\usage{
twRunDEMC(argsTwDEMCBatch, ..., argsReplaceZinit = list(), prevResRunCluster = NULL, restartFilename = NULL)
}
\arguments{
  \item{argsTwDEMCBatch}{Arguments passed to twDEMCBatch -> twDEMCInt -> fLogLik.
It is updated by \dots.
After update it must contain entries Zinit and fLogLik
It is further searched for entries nPops, and argsFLogLik, and resFLogLikX. The latter are initialized to defaults  \code{1,list(),character(0)} respectively if not found.   }
  \item{\dots}{further arguments passed to twDEMCBatch -> twDEMCInt -> fLogLik}
  \item{argsReplaceZinit}{further arguments passed to twDEMCBatch}
  \item{prevResRunCluster}{results of call to twRunDEMC, argument required to be called from runCluster.R}
  \item{restartFilename}{name of the file to store restart information, argument required to be called from runCluster.R}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





