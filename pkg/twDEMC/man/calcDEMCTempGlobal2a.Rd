\name{calcDEMCTempGlobal2a}
\alias{calcDEMCTempGlobal2a}
\title{calcDEMCTempGlobal2a}
\description{Calculating global temperature after the next batch for one population.}
\usage{
calcDEMCTempGlobal2a(resPop, diffLogDen, TLp, pAcceptTVar, iRun = getNGen(resPop), nGenBurnin, nRun, rHat0 = 1.08)
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
}

\value{list with components \itemize{
\item{TGlobal: numeric scalar: the global Temperature}
\item{nGenBurnin: recalculated burnin period}
}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





