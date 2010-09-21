\name{twCalcLogLikPar}
\alias{twCalcLogLikPar}
\title{twCalcLogLikPar}
\description{Invokes fLoglik with proposal in a parallel load balanced way.}
\usage{
twCalcLogLikPar(fLogLik, xProp, resFLogLikX = NULL, intResCompNames = character(0), argsFLogLik = list(), fLogLikScale = 1, debugSequential = FALSE, remoteDumpfileBasename = NULL, ...)
}
\arguments{
  \item{fLogLik}{the objective function}
  \item{xProp}{numeric matrix of proposals, columns: parameter vector components rows: cases}
  \item{resFLogLikX}{numeric matrix of result of fLogLik
colnames must contain intResCompNames 
rows: number of cases in xProp	}
  \item{intResCompNames}{character vector: names of results components of fLogLik that are used for internal Metropolis decisions }
  \item{argsFLogLik}{arguments passed to fLogLik}
  \item{fLogLikScale}{factor multiplied to the result of fLogLik}
  \item{debugSequential}{see \code{\link{sfFArgsApplyLB}}}
  \item{remoteDumpfileBasename}{see \code{\link{sfRemoteWrapper}}}
  \item{\dots}{further arguments passed to fLogLik}
}

\value{List with the following items \describe{
\item{logLik}{numeric vector: for each state: the sum of logLiks over all components, multiplied by fLogLikScale}
\item{resFLogLik}{numeric matrix: return components of fLogLik, one row for each state, columns: components }
}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twDEMCInt}}}

