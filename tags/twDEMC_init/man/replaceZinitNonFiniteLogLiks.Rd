\name{replaceZinitNonFiniteLogLiks}
\alias{replaceZinitNonFiniteLogLiks}
\title{replaceZinitNonFiniteLogLiks}
\description{Replaces states of Zinit that yield non-finite rLogLik by sampling other states.}
\usage{
replaceZinitNonFiniteLogLiks(Zinit, rLogLik)
}
\arguments{
  \item{Zinit}{initial states see InitDEMCzsp}
  \item{rLogLik}{tmp.rLogLik: calculated logLikelihoods for all the states in Zinit (rows cases and columns parameters)}
}
\details{In order for twDEMC to start up effectively, it is important that chains start from values, where the logLikelihood is finite}


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}
\code{\link{replaceZinitCases}}}

