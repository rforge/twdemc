\name{twDEMCSACont}
\alias{twDEMCSACont}
\title{twDEMCSACont}
\description{continuing simulated annealing DEMC based on previous result}
\usage{twDEMCSACont(mc, nGen = 512, nObs, ..., m0 = calcM0twDEMC(getNParms(mc), 
    getNChainsPop(mc)), controlTwDEMC = list(), debugSequential = FALSE, 
    restartFilename = NULL, ctrlBatch = list(nSamplesBatch = m0 * 
        8, pThinPast = 0.5), ctrlT = list(TFix = numeric(0), 
        TMax = numeric(0), TDecProp = 0.9), ctrlConvergence = list(maxRelTChangeCrit = 0.025, 
        maxLogDenDriftCrit = 0.3, gelmanCrit = 1.4, critSpecVarRatio = 20, 
        dumpfileBasename = NULL, maxThin = 0))}
\arguments{
  \item{mc}{result of twDEMCBlock}
  \item{nGen}{overall number of generations to add (within one batch provide ctrlBatch$nGenBatch)}
  \item{nObs}{integer vector (nResComp) specifying the number of observations for each result component}
  \item{\dots}{further argument to \code{\link{twDEMCBlockInt}}}
  \item{m0}{minimum number of samples in step for extending runs}
  \item{controlTwDEMC}{list argument to \code{\link{twDEMCBlockInt}} containing entry thin}
  \item{debugSequential}{set to TRUE to avoid parallel execution, good for debugging}
  \item{restartFilename}{filename to write intermediate results to}
  \item{ctrlBatch}{list of arguments controlling batch executions\describe{
\item{nSamplesBatch}{number of samples for one call to twDEMCStep (multiplied by thin to get generations)
<< default: set in a way that on average each population (assuming half are significant) is appended by 2*m0 samples}
\item{pThinPast}{in each batch thin the past to given fraction before appending new results}
}}
  \item{ctrlT}{list of arguments controlling Temperature decrease\describe{
\item{TFix}{numeric vector (nResComp) specifying a finite value for components with fixed Temperatue, and a non-finite for others}
\item{TMax}{numeric vector (nResComp) specifying a maximum temperature for result components.}
\item{TDecProp}{proportion of Temperature decrease: below one to diminish risk of decreasing Temperature too fast (below what is supported by other data streams)}
}}
  \item{ctrlConvergence}{list or arguments controlling check for convergence\describe{
\item{maxRelTChangeCrit}{if Temperature of the components changes less than specified value, the algorithm can finish}
\item{maxLogDenDriftCrit}{if difference between mean logDensity of first and fourth quartile of the sample is less than this value, we do not need further batches because of drift in logDensity}
\item{gelmanCrit}{do not change Temperature, if variance between chains is too high, i.e. Gelman Diag is above this value}
\item{critSpecVarRatio}{if proprotion of spectral Density to Variation is higher than this value, signal problems and resort to subspaces}
\item{dumpfileBasename}{scalar string: filename to dump stack before stopping. May set to "recover"}
\item{maxThin}{scalar positive integer: maximum thinning interval. If not 0 then thinning is increased on too high spectral density}
}}
}
\details{\describe{\item{diagnostics}{
in resulting twDEMC an entry diagnostics is returned. Its a list with diagnostics after each batch.
Its components have been calculated on the last halv of the chains.\describe{
\item{T0}{Temperature (variance inflation factor) scaled for one observation minus one}
\item{maxGelmanDiagPops}{maximum of within population gelman diagnostics}
\item{gelmanDiag}{gelman diag calculated between populations}
\item{maxSpecVarRatio}{ratio of spectral density to sample Variance (max over all parameters). This is a measure that increases with autocorrelation.
Last components represent the average logDensity of the last 8th of the chains.}
}

}}

\describe{\item{Adaptive thinning interval (\code{maxThin})}{ 
If spectral variance is larger than its critical ratio \code{ctrlConvergence$critSpecVarRatio}
then thinning interval needs to be increased.
If doubling the thinning interval is not larger than \code{ctrlConvergence$maxThin} then a new batch 
with higher thinning is attempted. Else the SADEMC quits with reporting an error with return values component \code{failureMsg}. 
}}}


\author{Thomas Wutzler}





