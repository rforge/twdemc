\name{twCalcLogDenPar}
\alias{twCalcLogDenPar}
\title{twCalcLogDenPar}
\description{Invokes fLogDen with proposal in a parallel load balanced way.}
\usage{twCalcLogDenPar(fLogDen, xProp, logDenCompX = NULL, intResComp = character(0), 
    argsFLogDen = list(), debugSequential = FALSE, remoteDumpfileBasename = NULL, 
    ...)}\arguments{
  \item{fLogDen}{the objective function}
  \item{xProp}{numeric matrix (nCases x nParm) of proposals}
  \item{logDenCompX}{all components of logDensity of xProp (result of fLogDen) 
colnames must contain intResCompNames 
rows: number of cases in xProp	}
  \item{intResComp}{character vector: names of results components of fLogDen that are used for internal Metropolis decisions }
  \item{argsFLogDen}{arguments passed to fLogDen}
  \item{debugSequential}{see \code{\link{sfFArgsApplyLB}}}
  \item{remoteDumpfileBasename}{see \code{\link{sfRemoteWrapper}}}
  \item{\dots}{further arguments passed to fLogDen}
}

\value{List with the following items \describe{
\item{logDen}{numeric vector: for each state: the sum of logDens over all components}
\item{logDenComp}{numeric matrix: return components of fLogDen, one row for each state, columns: components }
}}

\author{Thomas Wutzler}



\seealso{\code{\link{twDEMCBlockInt}}}
\examples{
data(twdemcEx1)
xProp <- stackChains(twdemcEx1)[,-1]

data(twLinreg1)
attach( twLinreg1 )
res <- twCalcLogDenPar(logDenGaussian,xProp
	,fModel=dummyTwDEMCModel		### the model function, which predicts the output based on theta 
	,obs=obs			### vector of data to compare with
	,invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
	thetaPrior = thetaTrue,	### the prior estimate of the parameters
	invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
	xval=xval
)
str(res)
}
