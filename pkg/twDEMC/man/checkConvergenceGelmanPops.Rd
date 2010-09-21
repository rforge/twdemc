\name{checkConvergenceGelmanPops}
\alias{checkConvergenceGelmanPops}
\title{checkConvergenceGelmanPops}
\description{Gelman RHat criterion applied to each population and between populations}
\usage{
checkConvergenceGelmanPops(res, burninFrac = 0.5, rHatMin = 1.1)
}
\arguments{
  \item{res}{see return value of \code{\link{twDEMCInt}} ($parms (d x nStep x nChain) )}
  \item{burninFrac}{fraction of the chain to be discarded (default 0.5)}
  \item{rHatMin}{rHat criterion, upper bound that is regarded as convergence}
}

\value{all rl <= criterion for each chain}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





