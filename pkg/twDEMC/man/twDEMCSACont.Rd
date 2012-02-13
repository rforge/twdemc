\name{twDEMCSACont}
\alias{twDEMCSACont}
\title{twDEMCSACont}
\description{continuing simulated annealing DEMC based on previous reslt}
\usage{twDEMCSACont(mc, nObs, nGen = 512, TFix = numeric(0), TMax = numeric(0), 
    nBatch = 4, maxRelTChange = 0.025, maxLogDenDrift = 0.3, 
    TDecProp = 0.9, gelmanCrit = 1.4, critSpecVarRatio = 20, 
    restartFilename = NULL, ...)}
\arguments{
  \item{mc}{result of twDEMCBlock}
  \item{nObs}{integer vector (nResComp) specifying the number of observations for each result component}
  \item{nGen}{number of generations in the initial batch, default 512}
  \item{TFix}{numeric vector (nResComp) specifying a finite value for components with fixed Temperatue, and a non-finite for others}
  \item{TMax}{numeric vector (nResComp) specifying a maximum temperature for result components.}
  \item{nBatch}{number of batches with recalculated Temperature}
  \item{maxRelTChange}{if Temperature of the components changes less than specified value, the algorithm can finish}
  \item{maxLogDenDrift}{if difference between mean logDensity of first and fourth quartile of the sample is less than this value, we do not need further batches because of drift in logDensity}
  \item{TDecProp}{proportion of Temperature decrease: below one to diminish risk of decreasing Temperature too fast (below what is supported by other data streams)}
  \item{gelmanCrit}{do not decrease Temperature, if variance between chains is too high, i.e. Gelman Diag is above this value}
  \item{critSpecVarRatio}{if proprotion of spectral Density to Variation is higher than this value, signal problems and resort to subspaces}
  \item{restartFilename}{filename to write intermediate results to}
  \item{\dots}{further argument to \code{\link{twDEMCBlockInt}}}
}



\author{Thomas Wutzler}





