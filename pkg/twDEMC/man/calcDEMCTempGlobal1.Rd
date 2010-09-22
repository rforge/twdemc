\name{calcDEMCTempGlobal1}
\alias{calcDEMCTempGlobal1}
\title{calcDEMCTempGlobal1}
\description{Calculating global temperature after the next batch for one population.}
\usage{
calcDEMCTempGlobal1(resPop, diffLogLik, TLp, pAcceptTVar, iRun = calcNGen(resPop), nGenBurnin, nRun, minPCompAcceptTempDecr = 0.16)
}
\arguments{
  \item{resPop}{twDEMC result (subChains of one population)}
  \item{diffLogLik}{numeric vector: Lp-La of the previous proposals}
  \item{TLp}{numeric scalar: max Temperature suggested by optimizing Lp  (from \code{\link{calcDEMCTempDiffLogLik3}}}
  \item{pAcceptTVar}{numeric scalar: Acceptance rate of temperatue dependent step (from \code{\link{calcDEMCTempDiffLogLik3}}}
  \item{iRun}{current generation: may be passed for efficiency}
  \item{nGenBurnin}{integer scalar: the number of Generations in burnin}
  \item{nRun}{integer scalar: the number of generations in next batch}
  \item{minPCompAcceptTempDecr}{
}
}
\details{This version either enforces Temp-Decrease complying to exponential decrease to 1 at nGenBurnin
or stays at temperature and prolonges nGenBurnin}
\value{list with components \itemize{
\item{TGlobal: numeric scalar: the global Temperature}
\item{nGenBurnin: recalculated burnin period}
}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





