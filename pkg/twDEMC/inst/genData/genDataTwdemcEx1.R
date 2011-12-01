#generating a typical twDEMC result
.nPop=2
.nChainPop=4
.thin=4
.nGenBurnin=30
.T0=10

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

Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChainPop=.nChainPop, nPop=.nPop)
dim(Zinit)

.nGen=100
#mtrace(twDEMCBlock.array)
#mtrace(twDEMCBlockInt)
twdemcEx1 <-  twDEMCBlock( Zinit, nGen=.nGen
	,dInfos=list(list(fLogDen=logDenGaussian, argsFLogDen=argsFLogDen))
	,nPop=.nPop
	,controlTwDEMC=list(thin=.thin)
	#,nGenBurnin=.nGenBurnin
	, T0=.T0
	#fLogDenScale=1		#default scale of logDenGaussian is already -1/2
	,doRecordProposals = TRUE
)
str(twdemcEx1)

.tmp.f <- function(){
	rescoda <- as.mcmc.list(twdemcEx1)
	plot(rescoda, smooth=FALSE)
}
save(twdemcEx1,file="data/twdemcEx1.RData")

detach( twLinreg1 )
