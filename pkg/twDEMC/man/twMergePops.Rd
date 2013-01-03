\name{twMergePops}
\alias{twMergePops}
\title{twMergePops}
\description{Merge subspaces that hold only few samples or have only low percentiles to neighboring populations}
\usage{twMergePops(mcApp, m0 = calcM0twDEMC(getNParms(mcApp), nChainPop), 
    minPSub = 0, pSubs = getNSamples(mcApp)/getNSamplesSpace(mcApp)[getSpacesPop(mcApp)], 
    nChainPop = getNChainsPop(mcApp), spaceInds = unique(getSpacesPop(mcApp)))}
\arguments{
  \item{mcApp}{twDEMCPop: holding populations to merge}
  \item{m0}{minimum number of samples}
  \item{minPSub}{minimal quantile of a sub-population}
  \item{pSubs}{numeric vector (nPops): percentiles of each population during last sampling (sum to 1)}
  \item{nChainPop}{integer scalar: number of chains per population, maybe passed for efficiency reasons}
  \item{spaceInds}{integer vector: set of space indices, maybe passed for efficiency reasons}
}

\value{list with entries
\item{mcApp}{twDEMCPop with merged populations}
\item{pSubs}{updated percentiles of the subPopulations. Sum to 1.}
\item{iPopsBefore}{integer vector (nPop): indices of the largest original population contributing to the new population}}

\author{Thomas Wutzler}





