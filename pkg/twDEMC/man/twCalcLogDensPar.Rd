\name{twCalcLogDensPar}
\alias{twCalcLogDensPar}
\title{twCalcLogDensPar}
\description{Invokes all fLogDens with proposal in a parallel way, taking care of intermediates between densities}
\usage{twCalcLogDensPar(dInfos,
 xProp, debugSequential = FALSE, remoteDumpfileBasename = NULL)}\arguments{
  \item{dInfos}{list describing the logDensities, see \code{\link{twDEMCBlockInt}}}
  \item{xProp}{numeric matrix (nCases x nParm) of proposals}
  \item{debugSequential}{see \code{\link{sfFArgsApplyLB}}}
  \item{remoteDumpfileBasename}{see \code{\link{sfRemoteWrapper}}}
}
\details{Does not take care of internal components: provides no argument logDenPrev}
\value{List with the following items \describe{
\item{logDen}{numeric vector: for each state: the sum of logDens over all components}
\item{logDenComp}{numeric matrix: return components of fLogDen, one row for each state, columns: components }
\item{logDenCompPos}{integer vector (nDenComp): index of the densitiy that provides the component }
}}

\author{Thomas Wutzler}



\seealso{\code{\link{twDEMCBlockInt}}, \code{\link{twCalcLogDenPar}}}
\examples{
# see test function twCalcLogDensPar in test case twCalcLogDenPar
}
