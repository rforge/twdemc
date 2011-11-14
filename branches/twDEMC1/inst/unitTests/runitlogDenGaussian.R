## twUtestF("logDenGaussian")

.setUp <-function () {
	.setUpDf <- within( list(),{
		fModel = dummyTwDEMCModel		# the model function, which predicts the output based on theta
		xval = 1:10		# argument needed by mofModelummy
		thetaTrue = c(a=2,b=5)	# the parameter vector
		sdy = xval^0.5
		obs = thetaTrue["a"] + thetaTrue["b"]*xval + rnorm(length(xval), sd=sdy)		### vector of data to compare with
	
		sdTheta= thetaTrue*0.05	# 5% error
		theta = thetaTrue + rnorm(length(thetaTrue),sd=sdTheta)
		
		#plot( xval, obs ); abline(2,5)
		invCovar = diag(1/sdy^2)		### the inverse of the Covariance of obs (its uncertainty)
	})
	attach(.setUpDf)
}

.tearDown <- function () {
	#detach(.setUpDf)
	detach()
}

test.noprior <- function (){
	#mtrace(logDenGaussian)
	res <- logDenGaussian( theta, fModel=dummyTwDEMCModel, obs=obs, invCovar=invCovar, xval=xval )
	msg <- as.character(res)
	checkTrue( is.numeric(res), msg )
	#checkEquals( length(res),1, msg)
	checkTrue( res["obs"] < 0, msg)
} 

test.prior <- function (){
	#thetaPrior = coef(lm(obs~xval))	### the prior estimate of the parameters
	#invCovarTheta = solve(summary(lm(obs~xval))$cov.unscaled)	### the inverse of the Covariance of the prior parameter estimates
	thetaPrior<- thetaTrue
	invCovarTheta <- diag(1/(sdTheta)^2)
	res <- logDenGaussian( theta, fModel=dummyTwDEMCModel, obs=obs, invCovar=invCovar, xval=xval, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta )
	msg <- as.character(res)
	checkTrue( is.numeric(res), msg )
	checkEquals( c("obs","parms"), names(res))
	checkTrue( all(res < 0))
} 

test.twoStepMetropolis <- function (){
	#thetaPrior = coef(lm(obs~xval))	### the prior estimate of the parameters
	#invCovarTheta = solve(summary(lm(obs~xval))$cov.unscaled)	### the inverse of the Covariance of the prior parameter estimates
	thetaPrior<- thetaTrue
	invCovarTheta <- diag(1/(sdTheta)^2)
	#mtrace(logDenGaussian)
	thetaProp=theta*1000	#some nearly unprobable combination
	res <- logDenGaussian( thetaProp, fModel=dummyTwDEMCModel, obs=obs, invCovar=invCovar, xval=xval, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta
		,logDenAccept=c(parms=-1e-10)	#provide a near one previous Density (near zero logDen)
	)
	msg <- as.character(res)
	checkTrue( is.numeric(res), msg )
	checkEquals( c(obs=NA, parms=-Inf), res )

	res <- logDenGaussian( thetaProp, fModel=dummyTwDEMCModel, obs=obs, invCovar=invCovar, xval=xval
		,logDenAccept=c(parms=-1e-10)	#provide a near one previous Density (near zero logDen)
	)
	msg <- as.character(res)
	checkTrue( is.numeric(res), msg )
	checkEquals( c("obs","parms"), names(res), msg)
	checkTrue( res["obs"] < 0, msg )
	checkTrue( res["parms"] == 0, msg)
	
} 

tmp.f <- function(){
	#library(debug)
	currentPackage="twDEMC"
	#mtrace(logDenGaussian)
	#twUtestF(logDenGaussian,"test.logDenGaussian")
	twUtestF(logDenGaussian)
	twUtestF()
}

