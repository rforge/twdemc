\name{divideTwDEMCStep}
\alias{divideTwDEMCStep}
\title{divideTwDEMCStep}
\description{run twDEMCBlock on subspaces}
\usage{divideTwDEMCStep(aTwDEMC, qPop = numeric(0), nGen = 512, controlTwDEMC = list(), 
    m0 = calcM0twDEMC(getNParms(aTwDEMC), getNChainsPop(aTwDEMC)), 
    nRowsMin, nRowsMax, minPSub, TEnd, ...)}
\arguments{
  \item{aTwDEMC}{former run of twDEMCBlockInt}
  \item{qPop}{numeric vector (nPop): probability of samples in subspace within entire space}
  \item{nGen}{the number of generations for twDEMCBlock}
  \item{controlTwDEMC}{list argument to \code{\link{twDEMCBlockInt}} containing entry thin}
  \item{m0}{minimum number of samples in step for extending runs}
  \item{nRowsMin}{number of rows in sub to sample at minimum, must contain enough samples for robust estimation of mean log-Density of subspace}
  \item{nRowsMax}{number of rows in sub to include as maximum mcNewMinN}
  \item{minPSub}{minimum proportion of a subsample in splitting}
  \item{TEnd}{numeric vector (nResComp) target temperature}
  \item{\dots}{further arguments to \code{\link{twDEMCBlockInt}}, such as TEnd}
}
\details{A first estimate of the proportion of samples from different subspaces
are the initial percentiles qp. 
The proportion of the samples from different subspaces is estimated
by the proportions of integrated density of the subspaces.
These proportions are esimated by average density over all samples multiplied by an
estimate proportional to the volume: the initial quantiles.}
\value{For each population, a list with entries
\item{resTwDEMC}{the runs with new samples of class twDEMCPops}
\item{pSubs}{the quantiles of samples within given subspace}
\item{subPercChange}{numeric vector: ratio of estimated proportion in limiting distribution to proportion of initial proportion}
\item{resTwDEMCMinN}{twDEMCPops with last (nRowsMin <- max( 2*m0, nSampleMin )) rows}
\item{nPopsApp}{numeric vectgor (nSpace): number of spaces for which samples are appended}}

\author{Thomas Wutzler}




\examples{
# see test function divideTwDEMCStep in test case divideTwDEMC 
}
