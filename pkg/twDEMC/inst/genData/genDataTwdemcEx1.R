#generating a typical twDEMC result
.nPops=2
.nChains=8
.thin=5
.nGenBurnin=30
.T0=10

data(twLinreg1)
attach( twLinreg1 )


argsFLogLik <- list(
	fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
	obs=obs,			### vector of data to compare with
	invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
	thetaPrior = thetaTrue,	### the prior estimate of the parameters
	invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
	xval=xval
)
do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))

Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=.nChains, nPops=.nPops)
dim(Zinit)

.nGen=100
twdemcEx1 <-  twDEMCBatch( Zinit, nGen=.nGen, 
	fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
	nPops=.nPops,
	controlTwDEMC=list(thin=.thin),
	nGenBurnin=.nGenBurnin, T0=.T0,
	fLogLikScale=1		#default scale logLikGaussian is -1/2
	,doRecordProposals = TRUE
)
str(twdemcEx1)

.tmp.f <- function(){
	rescoda <- as.mcmc.list(twdemcEx1)
	plot(rescoda)
}
save(twdemcEx1,file="data/twdemcEx1.RData")

detach( twLinreg1 )
