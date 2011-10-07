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
\examples{
data(twdemcEx1)
xProp <- stackChains(twdemcEx1$parms)

data(twLinreg1)
attach( twLinreg1 )
res <- twCalcLogLikPar(logLikGaussian,xProp
	,fModel=dummyTwDEMCModel		### the model function, which predicts the output based on theta 
	,obs=obs			### vector of data to compare with
	,invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
	thetaPrior = thetaTrue,	### the prior estimate of the parameters
	invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
	xval=xval
)
str(res)
resM <- matrix(res$logLik, ncol=dim(twdemcEx1$parms)[3])
str(resM)
}
