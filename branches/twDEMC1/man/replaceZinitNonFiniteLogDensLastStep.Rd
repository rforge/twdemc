\name{replaceZinitNonFiniteLogDensLastStep}
\alias{replaceZinitNonFiniteLogDensLastStep}
\title{replaceZinitNonFiniteLogDensLastStep}
\description{Replaces states of last step, i.e. column of Zinit that yield non-finite rLogDen by sampling other states.}
\usage{
replaceZinitNonFiniteLogDensLastStep(Zinit, fLogDen, nPops = 1, iStep = dim(Zinit)[2], maxSteps = 16, ...)
}
\arguments{
  \item{Zinit}{initial states see InitDEMCzsp}
  \item{fLogDen}{the logDen Function}
  \item{nPops}{number of populations. States are only choosen from same population}
  \item{iStep}{the step for which to replace nonfinite yielding parameters.}
  \item{maxSteps}{
}
  \item{\dots}{arguments to \code{\link{twCalcLogDenPar}}}
}

\value{list with components \describe{
\item{Zinit}{initial states see InitDEMCzsp}
\item{rLogDen, rLogDenComp, logDenComp}{numeric matrix: results of \code{\link{twCalcLogDenPar}} for last row for all chains}
}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}}

