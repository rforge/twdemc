\name{replaceZinitNonFiniteLogDens}
\alias{replaceZinitNonFiniteLogDens}
\title{replaceZinitNonFiniteLogDens}
\description{Replaces states of Zinit that yield non-finite rLogDen by sampling other states.}
\usage{
replaceZinitNonFiniteLogDens(Zinit,
logDen)
}
\arguments{
  \item{Zinit}{initial states see InitDEMCzsp}
  \item{logDen}{numeric matrix (nCases x nChains): calculated logDensitys for all the states in Zinit. If it is a vector then it is reshaped.}
}
\details{In order for twDEMC to start up effectively, it is important that chains start from values, where the logDensity is finite}
\value{Zinit, with several cols (parameter vectors) replaced by other cols}

\author{Thomas Wutzler}



\seealso{\code{\link{initZtwDEMCNormal}}
\code{\link{replaceZinitCases}}}
\examples{
data(twdemcEx1)
data(twLinreg1)
attach( twLinreg1 )
#mtrace(initZtwDEMCNormal)
Zinit <- initZtwDEMCNormal( theta0, diag(4*sdTheta^2), nChainPop=8, nPop=2)
dim(Zinit)
res <- res0 <- twCalcLogDenPar(logDenGaussian, stackChains(Zinit) # chains stack to calculate in parallel
	,fModel=dummyTwDEMCModel		### the model function, which predicts the output based on theta 
	,obs=obs				### vector of data to compare with
	,invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
	thetaPrior = thetaTrue,	### the prior estimate of the parameters
	invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
	xval=xval
)$logDen
plot(density(res))
res[res < -30] <- NA
resM <- matrix(res, ncol=dim(Zinit)[3] ) # restacking to chains work because it is the last dimension
set.seed(0815)
Zinit2 <- replaceZinitNonFiniteLogDens(Zinit, resM)
set.seed(0815)
Zinit3 <- replaceZinitNonFiniteLogDens(Zinit, res)
identical(Zinit2,Zinit3)

res2 <- twCalcLogDenPar(logDenGaussian, stackChains(Zinit2) # chains stack to calculate in parallel
	,fModel=dummyTwDEMCModel		### the model function, which predicts the output based on theta 
	,obs=obs			### vector of data to compare with
	,invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
	thetaPrior = thetaTrue,	### the prior estimate of the parameters
	invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
	xval=xval
)$logDen
d2 <- density(res2)
lines( d2$x, d2$y * max(density(res0)$y)/max(d2$y), col="blue")
}
