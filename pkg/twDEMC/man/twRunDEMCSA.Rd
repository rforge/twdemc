\name{twRunDEMCSA}
\alias{twRunDEMCSA}
\title{twRunDEMCSA}
\description{wrapper to twDEMCSA compliant to runCluster.R}
\usage{twRunDEMCSA(argsTwDEMCSA, ..., prevResRunCluster = NULL, restartFilename = NULL)}
\arguments{
  \item{argsTwDEMCSA}{Arguments passed to twDEMCSA -> twDEMCBlockInt
It is updated by \dots.}
  \item{\dots}{further arguments passed to twDEMCSA -> twDEMCBlockInt}
  \item{prevResRunCluster}{results of call to twRunDEMCSA, argument required to be called from runCluster.R}
  \item{restartFilename}{name of the file to store restart information, argument required to be called from runCluster.R}
}



\author{Thomas Wutzler}





