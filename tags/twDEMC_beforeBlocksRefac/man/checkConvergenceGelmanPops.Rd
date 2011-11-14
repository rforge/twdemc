\name{checkConvergenceGelmanPops}
\alias{checkConvergenceGelmanPops}
\title{checkConvergenceGelmanPops}
\description{Gelman RHat criterion applied to each population and between populations}
\usage{
checkConvergenceGelmanPops(aTwDEMC, burninFrac = aTwDEMC$nGenBurnin/getNGen(aTwDEMC), rHatMin = 1.1)
}
\arguments{
  \item{aTwDEMC}{see return value of \code{\link{twDEMCInt}} ($parms (d x nStep x nChain) )}
  \item{burninFrac}{fraction of the chain to be discarded (default 0.5)}
  \item{rHatMin}{rHat criterion, upper bound that is regarded as convergence}
}

\value{all rl <= criterion for each chain}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(twdemcEx1)
aTwDEMC <- twdemcEx1
checkConvergenceGelmanPops(twdemcEx1)
}
