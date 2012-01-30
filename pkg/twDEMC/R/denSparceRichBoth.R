denSparce <- function(
	### Example of using two different logDensity functions: density of sparce observations
	theta		
	,twTwoDenEx=twTwoDenEx1
	,theta0=twTwoDenEx$thetaTrue
	,...
){
	theta0[names(theta)] <- theta
	pred <- twTwoDenEx$fModel(theta0, xSparce=twTwoDenEx$xSparce, xRich=twTwoDenEx$xRich, ...)
	misfit <- twTwoDenEx$obs$y1 - pred$y1
	-1/2 * sum((misfit/twTwoDenEx$sdObs$y1)^2)
}

denSparcePrior <- function(
	### Same as denSpace but returning two components: prior and misfit
	theta		
	,twTwoDenEx=twTwoDenEx1
	,theta0=twTwoDenEx$thetaTrue
	,thetaPrior = NULL		##<< the prior estimate of the parameters
	,invCovarTheta = NULL	##<< the inverse of the Covariance of the prior parameter estimates
	,...
){
	theta0[names(theta)] <- theta
	pred <- twTwoDenEx$fModel(theta0, xSparce=twTwoDenEx$xSparce, xRich=twTwoDenEx$xRich, ...)
	misfit <- twTwoDenEx$obs$y1 - pred$y1
	logDenPropParms <- if( !is.null(thetaPrior) ){
			tmp.diffParms <- theta - thetaPrior
			if( is.matrix(invCovarTheta) )
				t(tmp.diffParms) %*% invCovarTheta %*% tmp.diffParms
			else
				sum(tmp.diffParms^2 / invCovarTheta)	# assume invCovar to be independent variances
		} else 0
	-1/2 * c( parmsSparce=logDenPropParms , obsSparce=  sum((misfit/twTwoDenEx$sdObs$y1)^2) )
}


denRich <- function(
	### Example of using two different logDensity functions: density of data-rich observations
	theta
	,twTwoDenEx=twTwoDenEx1
	,theta0=twTwoDenEx$thetaTrue
	,...
){
	theta0[names(theta)] <- theta
	pred <- twTwoDenEx$fModel(theta0, xSparce=twTwoDenEx$xSparce, xRich=twTwoDenEx$xRich,...)
	misfit <- twTwoDenEx$obs$y2 - pred$y2
	-1/2 * sum((misfit/twTwoDenEx$sdObs$y2)^2)
}

denBoth <- function(
	### Example of using two different logDensity functions: comparison of combining denSparce and denRich into one function.
	theta
	,twTwoDenEx=twTwoDenEx1
	,weights=c(1,1,1)			# weights for the two data streams
	,theta0=twTwoDenEx$thetaTrue
	,thetaPrior = NULL		##<< the prior estimate of the parameters
	,invCovarTheta = NULL	##<< the inverse of the Covariance of the prior parameter estimates
	,...
){
	logDenPropParms <- if( !is.null(thetaPrior) ){
			tmp.diffParms <- theta - thetaPrior
			if( is.matrix(invCovarTheta) )
				t(tmp.diffParms) %*% invCovarTheta %*% tmp.diffParms
			else
				sum(tmp.diffParms^2 / invCovarTheta)	# assume invCovar to be independent variances
		} else 0
	theta0[ names(theta) ] <- theta
	pred <- twTwoDenEx$fModel(theta0, xSparce=twTwoDenEx$xSparce, xRich=twTwoDenEx$xRich, ...)
	misfit <- twTwoDenEx$obs$y1 - pred$y1
	d1 <- sum((misfit/twTwoDenEx$sdObs$y1)^2)
	misfit <- twTwoDenEx$obs$y2 - pred$y2
	d2 <- sum((misfit/twTwoDenEx$sdObs$y2)^2)
	db <- -1/2 * c(parms=logDenPropParms, y1=d1,y2=d2)*weights
	db
}
attr(denBoth,"ex") <- function(){
	data(twTwoDenEx1)
	denBoth(twTwoDenEx1$thetaTrue)
	denBoth(twTwoDenEx1$thetaTrue, thresholdCovar = 0.3)
}
