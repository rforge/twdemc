data(twLinreg1)

denSigmaEx1Obs <- function(
	theta				 ##<< named numeric vector a,b,logSigma2
	,twLinreg=twLinreg1 ##<< list with components xval and obs 
){
	pred <- twLinreg$fModel(theta[1:2],twLinreg$xval)
	resid <- pred - twLinreg$obs
	varObs <- exp(theta[3])
	-1/2 * sum(resid^2)/varObs
}
attr(denSigmaEx1Obs,"ex") <- function(){
	data(twLinreg1)
	#mtrace(denSigmaEx1Obs)
	logSigma2 <- 2*log( mean(twLinreg1$sdObs) )
	denSigmaEx1Obs( c( twLinreg1$theta0, logSigma2=logSigma2 ) )
	
	#likelihood profile of a and b conditional on given logSigma2
	lmDummy <- with(twLinreg1, lm( obs ~ xval, weights=1/sdObs^2))		# results without priors
	(.expTheta <- coef(lmDummy))
	(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
	(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=c("a","b")) )
	nGrid <- 50
	aGrid <- .expTheta[1] + .expSdTheta[1]*seq(-2, 2, length.out=nGrid)	
	bGrid <- .expTheta[2] + .expSdTheta[2]*seq(-2, 2, length.out=nGrid)
	abGrid <- expand.grid(a=aGrid,b=bGrid)
	Zinit <- cbind(abGrid, logSigma2=logSigma2)
	logLikab <- apply(Zinit, 1, denSigmaEx1Obs  )
	logLikabm <- matrix(logLikab, nrow=nGrid)
	logLikabm[ logLikabm < max(logLikabm)-3 ] <- NA	# do not color regions with very low likelihood
	image( aGrid, bGrid, logLikabm, col=rev(heat.colors(80)) )
}


denSigmaEx1Sigma <- function(
	theta				 ##<< named numeric vector a,b,logSigma2
	,twLinreg=twLinreg1 ##<< list with components fModel, xval and obs
	,varSigma2=var( twLinreg1$sdObs^2 )	##<< variance of the normal distribution of uncertainty of sigma2 
){
	pred <- twLinreg$fModel(theta[1:2],twLinreg$xval)
	resid <- twLinreg$obs - pred
	sigma2Obs <- var(resid)
	-1/2 * (exp(theta[3]) - sigma2Obs)^2/varSigma2 
}
attr(denSigmaEx1Sigma,"ex") <- function(){
	data(twLinreg1)
	#mtrace(denSigmaEx1Sigma)
	logSigma2 <- 2*log( mean(twLinreg1$sdObs) )
	denSigmaEx1Sigma( c( twLinreg1$theta0, logSigma2=logSigma2 ) )
	
	#likelihood profile of varSigma conditional on the true parameters
	Zinit <- cbind( matrix(twLinreg1$thetaTrue,nrow=41,ncol=2,byrow=TRUE, dimnames=list(case=NULL,par=names(twLinreg1$thetaTrue)))
	, logSigm2=logSigma2*seq(0.8,1.2,length.out=41))
	head(Zinit)
	logLikSigma <- apply( Zinit, 1, denSigmaEx1Sigma, varSigma2=var( twLinreg1$sdObs^2 ) )
	bo <- logLikSigma > max(logLikSigma)-3
	plot( logLikSigma[bo] ~ exp(Zinit[bo,3]/2) )
	abline(h=-2)
}


logSigma2 <- 2*log( mean(twLinreg1$sdObs) )		# expected sigma2
varLogSigma2=1#var( 2*log(twLinreg1$sdObs) )		# variance for initialization

lmDummy <- with(twLinreg1, lm( obs ~ xval, weights=1/sdObs^2))		# results without priors
(.expTheta <- coef(lmDummy))
(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=c("a","b")) )


.nPop=2
.nChainPop=4
ZinitPops <- with(twLinreg1, initZtwDEMCNormal( c(theta0,logSigma2), diag(c(varTheta,varLogSigma2)), nChainPop=.nChainPop, nPop=.nPop))
#dim(ZinitPops)
#head(ZinitPops[,,1])
pops <- pops0 <- list(
	pop1 <- list(
		Zinit = ZinitPops[1:3,,1:4,drop=FALSE]	# the first population with less initial conditions
		,nGen=8
	),
	pop2 <- list(
		Zinit = ZinitPops[,,5:8,drop=FALSE]	# the first population with less initial conditions
		,nGen=100
		,T0=10
	)
)
#tmp <- .checkPop(pops[[1]])
# both fLogDen compare the same model against the same observations but use different priors
dInfoDefault <- list(
	fLogDen=logDenGaussian
)
dInfos0 <- dInfox <- list(
	rich = within(dInfoDefault,{
			argsFLogDen = list(
				fModel=rich$fModel,		### the model function, which predicts the output based on theta 
				obs=rich$obs,			### vector of data to compare with
				invCovar=rich$varObs,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
				xval=rich$xval
			)
		})
	,sparce = within(dInfoDefault,{
			argsFLogDen = list(
				fModel=sparce$fModel,		### the model function, which predicts the output based on theta 
				obs=sparce$obs,			### vector of data to compare with
				invCovar=sparce$varObs,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
				xval=sparce$xval
			)
		})
	,both = list( fLogDen=denBoth, argsFLogDen = list(twTwoDenEx1=twTwoDenEx1) )
)
#mtrace(logDenGaussian)
do.call( logDenGaussian, c(list(theta=rich$theta0), dInfos0$rich$argsFLogDen))
do.call( logDenGaussian, c(list(theta=sparce$theta0), dInfos0$sparce$argsFLogDen))
do.call( denBoth, c(list(theta=rich$theta0), dInfos0$both$argsFLogDen))


#fit only to rich data
blocks <- list(
	rich <- list(compPos=c("a","b"), dInfoPos="rich")
)
.nGen=100
.thin=4
#mtrace(.updateBlockTwDEMC)
#mtrace(.updateBlocksTwDEMC)
#mtrace(.updateIntervalTwDEMCPar)
#mtrace(twDEMCBlockInt)
resRichOnly <- res <- twDEMCBlock( pops=pops, dInfos=dInfos0, blocks=blocks, nGen=.nGen, controlTwDEMC=list(thin=.thin) )
str(res$pops[[2]])















detach( twTwoDenEx1 )