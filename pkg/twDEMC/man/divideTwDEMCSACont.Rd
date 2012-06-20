\name{divideTwDEMCSACont}
\alias{divideTwDEMCSACont}
\title{divideTwDEMCSACont}
\description{run twDEMCBlock on subspaces with decreasing Temperature}
\usage{
divideTwDEMCSACont(mc,
nGen = 512,
nObs,
iPopsDoSplit = integer(0),
...,
controlTwDEMC = list(),
doRecordProposals = FALSE,
m0 = calcM0twDEMC(getNParms(mc),
getNChainsPop(mc)),
debugSequential = FALSE,
restartFilename = NULL,
ctrlBatch = list(nGenBatch = m0 * (1/2/ctrlSubspaces$minPSub) * controlTwDEMC$thin * 2,
nSampleMin = 32,
pThinPast = 0.5),
ctrlT = list(TFix = numeric(0),
TMax = numeric(0),
TDecProp = 0.9),
ctrlConvergence = list(maxSubPercChangeCrit = 1.6,
maxRelTChangeCrit = 0.025,
maxLogDenDriftCrit = 0.3,
gelmanCrit = 1.4,
critSpecVarRatio = 10,
pCheckSkipPart = 0.5,
minNSampleCheck = 25,
dumpfileBasename = NULL),
ctrlSubspaces = list(minPSub = 0.1,
maxNSample = max(256,
2 * m0 * getNChainsPop(mc)),
argsFSplit = list()))
}
\arguments{
  \item{mc}{former run of twDEMCBlockInt}
  \item{nGen}{the number of generations for twDEMCBlock}
  \item{nObs}{integer vector (nResComp) specifying the number of observations for each result component}
  \item{iPopsDoSplit}{populations given here will definitely be splitted (no matter of other criteria)}
  \item{\dots}{arguments passed to \code{\link{twDEMCBlockInt}}}
  \item{controlTwDEMC}{list argument to \code{\link{twDEMCBlockInt}} containing entry thin}
  \item{doRecordProposals}{if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y}
  \item{m0}{minimum number of samples in step for extending runs

minimum number of samples in step for extending runs}
  \item{debugSequential}{set to TRUE to avoid parallel execution, good for debugging}
  \item{restartFilename}{filename to write intermediate results to}
  \item{ctrlBatch}{list of arguments controlling batch executions\describe{
\item{nGenBatch}{number of generations for one call to twDEMCStep
<< default: set in a way that on average each population (assuming half are significant) is appended by 2*m0 samples}
\item{nSampleMin}{minimum number of samples in each population within batch so that calculation of average density is stable}
\item{pThinPast}{in each batch thin the past to given fraction before appending new results}
}}
  \item{ctrlT}{list of arguments controlling Temperature decrease\describe{
\item{TFix}{numeric vector (nResComp) specifying a finite value for components with fixed Temperatue, and a non-finite for others}
\item{TMax}{numeric vector (nResComp) specifying a maximum temperature for result components.}
\item{TDecProp}{proportion of Temperature decrease: below one to diminish risk of decreasing Temperature too fast (below what is supported by other data streams)}
}}
  \item{ctrlConvergence}{list or arguments controlling check for convergence\describe{
\item{maxSubPercChangeCrit}{if all subPercChange of all sub-populations are below this value in two last batches, may assume convergence and skip further batches}
\item{maxRelTChangeCrit}{if Temperature of the components changes less than specified value, the algorithm can finish}
\item{maxLogDenDriftCrit}{if difference between mean logDensity of first and fourth quartile of the sample is less than this value, we do not need further batches because of drift in logDensity}
\item{gelmanCrit}{do not change Temperature, if variance between chains is too high, i.e. Gelman Diag is above this value}
\item{critSpecVarRatio}{if proprotion of spectral Density to Variation is higher than this value, split the space}
\item{pCheckSkipPart}{when checking each population for convergence problems, skip this proportion (kind of burnin)}
\item{minNSampleCheck}{if a population has less samples within a chains, skip the check because it is too uncertain}
\item{dumpfileBasename}{scalar string: filename to dump stack before stopping. May set to "recover"}
}}
  \item{ctrlSubspaces}{list of arguments controlling splitting and merging of subspaces\describe{
\item{minPSub}{minimal proportion of a sub-population}
\item{maxNSample}{if given a value, then looking for subspaces is done on a subsample of given length for efficiency (see \code{\link{getSubSpaces}})}
\item{argsFSplit}{further arguments to \code{\link{findSplit}}}
}}
}

\value{twDEMC with additional list entries
\item{nGen}{scalar integer: number of completed generations}
\item{pSubsBatch}{numeric vector (nPop): proportion of the populations within space during last batch}
\item{subPercChange}{numeric vector(nPop): relative change of proportions during last batch}
\item{relTChange}{numeric vector (nDen): calculated relative change of calculated new temperature for the next batch}
\item{specVarRatioPop}{numeric vector (nPop): ratio of spectral density to variance}
\item{TGlobal}{Temperature of non-fixed components}
\item{args}{calling arguments to provide restart capability}}

\author{Thomas Wutzler}





