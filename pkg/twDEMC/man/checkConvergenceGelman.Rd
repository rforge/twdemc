\name{checkConvergenceGelman}
\alias{checkConvergenceGelman}
\title{checkConvergenceGelman}
\description{Gelman RHat criterion applied to twDEMC result assuming all chains one population.}
\usage{
checkConvergenceGelman(res, burninFrac = 0.5, rHatMin = 1.1)
}
\arguments{
  \item{res}{see return value of \code{\link{twDEMCInt}} ($parms (d x nStep x nChain) )}
  \item{burninFrac}{fraction of the chain to be discarded (default 0.5)}
  \item{rHatMin}{rHat criterion, upper bound that is regarded as convergence}
}
\details{see Gelman04 (twutz:Gelman04_3#Inference_and_assessing_convergence)}
\value{all rl <= criterion for each chain}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





