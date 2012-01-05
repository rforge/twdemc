\name{replaceZinitNonFiniteLogDensLastStep}
\alias{replaceZinitNonFiniteLogDensLastStep}
\title{replaceZinitNonFiniteLogDensLastStep}
\description{Replaces states of last step, i.e. column of Zinit that yield non-finite rLogDen by sampling other states.}
\usage{
replaceZinitNonFiniteLogDensLastStep(Zinit,
fLogDen,
nPop = 1,
iStep = nrow(Zinit),
maxSteps = 16,
...)
}
\arguments{
  \item{Zinit}{initial states see InitDEMCzsp}
  \item{fLogDen}{the logDen Function}
  \item{nPop}{number of populations. States are only choosen from same population}
  \item{iStep}{the step for which to replace nonfinite yielding parameters.}
  \item{maxSteps}{maximum number of iterations (when last row was replaced by another row yielding non-finite loglik)}
  \item{\dots}{arguments to \code{\link{twCalcLogDenPar}}}
}

\value{list with components
\item{Zinit}{argument \code{Zinit} with some states in last row replaced by other states from Zinit.}
\item{rLogDen}{numeric vector: results of \code{\link{twCalcLogDenPar}} for last row for all chains}
\item{logDenComp}{numeric matrix: results of \code{\link{twCalcLogDenPar}} for last row for all chains
}}

\author{Thomas Wutzler}



\seealso{\code{\link{initZtwDEMCNormal}}}

