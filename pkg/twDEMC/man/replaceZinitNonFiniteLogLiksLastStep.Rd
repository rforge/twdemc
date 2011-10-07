\name{replaceZinitNonFiniteLogLiksLastStep}
\alias{replaceZinitNonFiniteLogLiksLastStep}
\title{replaceZinitNonFiniteLogLiksLastStep}
\description{Replaces states of last step, i.e. column of Zinit that yield non-finite rLogLik by sampling other states.}
\usage{
replaceZinitNonFiniteLogLiksLastStep(Zinit, fLogLik, nPops = 1, iStep = dim(Zinit)[2], maxSteps = 16, ...)
}
\arguments{
  \item{Zinit}{initial states see InitDEMCzsp}
  \item{fLogLik}{the logLik Function}
  \item{nPops}{number of populations. States are only choosen from same population}
  \item{iStep}{the step for which to replace nonfinite yielding parameters.}
  \item{maxSteps}{
}
  \item{\dots}{arguments to \code{\link{twCalcLogLikPar}}}
}

\value{list with components \describe{
\item{Zinit}{initial states see InitDEMCzsp}
\item{rLogLik, rLogLikComp, resFLogLik}{numeric matrix: results of \code{\link{twCalcLogLikPar}} for last row for all chains}
}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}}

