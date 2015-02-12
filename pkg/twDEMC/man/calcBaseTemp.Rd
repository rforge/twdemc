\name{calcBaseTemp}
\alias{calcBaseTemp}
\title{calcBaseTemp}
\description{calculate Base temperature from given temperatures of data streams and corresponding number of observations.}
\usage{calcBaseTemp(temp, nObs, TFix = rep(NA_real_, length(temp)), 
    iNonFixTemp = which(!is.finite(TFix)))}\arguments{
  \item{temp}{numeric vector (nResultComp): temperature, i.e. variance inflation factor}
  \item{nObs}{integer vector (nResultComp): number of observations}
  \item{TFix}{numeric vector (nResultComp): fixed temperature for components, non-finite for those with varying temperature
, alternatively a named vector listing only the components with fixed temperatures (temp must have names then too)}
  \item{iNonFixTemp}{integer vector: index of result components, which Temperature is fixed}
}
\details{The Temperature, i.e. variance inflation factor, scales with the number of observations by
T = 1 + (T0 - 1) * n.  
A common base temperature is the maximium across calculated stream base temperatures}
\value{numeric scalar: The temperature, i.e. variance inflation factor, for a single observation.}

\author{Thomas Wutzler}



\seealso{\link{calcStreamTemp}, \link{twDEMC}

\code{\link{calcTemperatedLogDen.default}}}
\examples{
data(twdemcEx1)   
.nObs <- c(parms=getNParms(twdemcEx1), obs=length(twdemcEx1$dInfos[[1]]$argsFLogDen$obs) )
.T <- getCurrentTemp(twdemcEx1)
calcBaseTemp( .T, .nObs[names(.T)], TFix=c(parms=1) )
}
