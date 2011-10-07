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
\value{Zinit, with several cols (parameter vectors) replaced by other cols}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}
\code{\link{replaceZinitCases}}}
\examples{
data(twLinreg1)
attach( twLinreg1 )
mtrace(initZtwDEMCNormal)
Zinit <- initZtwDEMCNormal( theta0, diag(4*sdTheta^2), nChains=8, nPops=2)
dim(Zinit)
res <- twCalcLogLikPar(logLikGaussian, stackChains(Zinit)
	,fModel=dummyTwDEMCModel		### the model function, which predicts the output based on theta 
	,obs=obs			### vector of data to compare with
	,invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
	thetaPrior = thetaTrue,	### the prior estimate of the parameters
	invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
	xval=xval
)$logLik
plot(density(res))
res[res< -30] <- NA
resM <- matrix(res, ncol=dim(twdemcEx1$parms)[3])
Zinit2 <- replaceZinitNonFiniteLogLiks(Zinit, resM)
}
