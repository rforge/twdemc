\name{twDEMCBlock.array}
\alias{twDEMCBlock.array}
\title{twDEMCBlock array}
\description{Initialize \code{\link{twDEMCBlockInt}} by array of initial population and remove those generations from results afterwards}
\usage{
\method{twDEMCBlock}{array}(x,
...,
nPop = 1,
X = NULL,
logDenCompX = NULL,
upperParBounds = vector("list",
nPop),
lowerParBounds = vector("list",
nPop))
}
\arguments{
  \item{x}{initial population: a numeric array (M0 x d x nChain) see details in \code{\link{twDEMCBlockInt}}}
  \item{\dots}{further arguments to \code{\link{twDEMCBlockInt}}}
  \item{nPop}{number of populations in x}
  \item{X}{initial state (nParm x nChain)}
  \item{logDenCompX}{numeric matrix (nResComp x nChains) initial state}
  \item{upperParBounds}{list of named numeric vectors: giving upper parameter bounds for each population 
for exploring subspaces of the limiting distribution, see details
\cr Alternatively a single numeric vector can be supplied, which is replicated for each population.}
  \item{lowerParBounds}{similar to upperParBounds}
}



\author{Thomas Wutzler}




\examples{
data(twLinreg1)
attach( twLinreg1 )

argsFLogDen <- list(
	fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
	obs=obs,			### vector of data to compare with
	invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
	thetaPrior = thetaTrue,	### the prior estimate of the parameters
	invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
	xval=xval
)
do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))

.nPop = 2
Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChainPop=4, nPop=.nPop)
dim(Zinit)

.nGen=100
#nGen=3
#mtrace(twDEMC.array)
#mtrace(.updateIntervalTwDEMCPar)
#mtrace(twDEMCBlockInt)
tmp1 <- tmp <-  twDEMCBlock( Zinit, nPop=.nPop
	,dInfos=list(list(fLogDen=logDenGaussian, argsFLogDen=argsFLogDen))
	,nGen=.nGen 
)
plot( as.mcmc.list(tmp), smooth=FALSE )
tmp2 <- tmp <- twDEMCBlock( tmp1, nGen=200 )

str(tmp)


}
