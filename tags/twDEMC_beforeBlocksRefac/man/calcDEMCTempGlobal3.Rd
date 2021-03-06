\name{calcDEMCTempGlobal3}
\alias{calcDEMCTempGlobal3}
\title{calcDEMCTempGlobal3}
\description{Calculating global temperature and adjusted burnin period after the next batch for one population.}
\usage{
calcDEMCTempGlobal3(resPop, diffLogDen, TLp, pAcceptTVar, iRun = getNGen(resPop), nGenBurnin, nRun, rHat0 = 1.2, nSampleGelmanDiag = 64, nSampleTrend = 512)
}
\arguments{
  \item{resPop}{twDEMC result (subChains of one population)}
  \item{diffLogDen}{numeric vector: Lp-La of the previous proposals}
  \item{TLp}{numeric scalar: max Temperature suggested by optimizing Lp  (from \code{\link{calcDEMCTempDiffLogDen3}}}
  \item{pAcceptTVar}{numeric scalar: Acceptance rate of temperatue dependent step (from \code{\link{calcDEMCTempDiffLogDen3}}}
  \item{iRun}{current generation: may be passed for efficiency}
  \item{nGenBurnin}{integer scalar: the number of Generations in burnin}
  \item{nRun}{integer scalar: the number of generations in next batch}
  \item{rHat0}{rHat value for which to not change the burnin period}
  \item{nSampleGelmanDiag}{sample size of the appropriately thinned sample that will be detrended to check Gelman diagnostics before decreasing temperature.
When working in high dimensional parameter space, you may encounter false high gelman diagnostics. Then you can try increasing this parameter,
but be aware of increasing computational demand.}
  \item{nSampleTrend}{size of the end of the appropriately thinned logDen sample that will be checked for a trend.}
}

\value{list with components \itemize{
\item{TGlobal: numeric scalar: the global Temperature}
\item{nGenBurnin: recalculated burnin period}
}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





