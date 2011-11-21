logDenGaussian <- function(
	### Invokes the model and calculate an unnormalized logDensity (-1/2*misfit) assuming (multivariate) Gaussian errors in both data in priors. 
	theta,			##<< the parameter vector to be optimized (may be more than updated, if used with blocks).
	logDenAccept=numeric(0),	##<< scalar: logDen for parms from revious run for two step Metropolis decision
	metropolisStepTemp=c(parms=1),		##<< numeric named vector: the temperature for internal metropolis step
	..., 			##<< any other arguments passed to fModel
	fModel,			##<< the model function, which predicts the output based on theta 
	theta0=theta,	##<< parameter vector, first argument to fModel. Before invocation components theta overwrite theta0 
	obs,			##<< vector of data to compare with
	invCovar,		##<< the inverse of the Covariance of obs (its uncertainty)
		##<< alternatively a vector of variances (diagonal covariance matrix) can be supplied and calculation is much more efficient
	thetaPrior = NULL,	##<< the prior estimate of the parameters
	invCovarTheta = NULL,	##<< the inverse of the Covariance of the prior parameter estimates
		##<< alternatively a vector of variances (diagonal covariance matrix) can be supplied and calculation is much more efficient
	namesTheta=NULL, ##<< names assigned to theta (if not NULL), before invoking mofModel
	blockIndices=NULL,	##<< integer vector: index of the components in theta and theta0 that should be regarded in this block 
	scale=-1/2 	 		##<< factor to mulitply the misfit (e.g. -1/2 to obtain the unnormalized logDensity)
){
	# logDenGaussian
	##seealso<<   
	## \code{\link{twDEMCInt}}
	## \code{\link{dummyTwDEMCModel}}
	if( !is.null(namesTheta ))
		names(theta) <- namesTheta	##details<<
	theta0[names(theta)] <- theta
	
	# if blockIndices are given, only evaluate prior of subset of parameter components	
	thetaBlock <- if( 0==length(blockIndices) ) theta else theta[blockIndices]
	
	##details<< 
	## If thetaPrior is not specified (NULL) then no penalty is assigned to parameters.
	logDenPropParms <- if( !is.null(thetaPrior) ){
			tmp.diffParms <- thetaBlock - thetaPrior
			if( is.matrix(invCovarTheta) )
				t(tmp.diffParms) %*% invCovarTheta %*% tmp.diffParms
			else
				sum(tmp.diffParms^2 / invCovarTheta)	# assume invCovar to be independent variances
		} else 0
	##details<<
	## Supports a two-step Metropolis descision. If \code{logDenAccept["parms"]} is provided, 
	## then a Metropolis descision is done based only on the parameters.
	## If it fails, then \code{c(obs=NA, parms=-Inf)} is returned. 
	## The possible costly evaluation of fModel is avoided.
	if( !is.na((logDenXParms <- logDenAccept["parms"])) & (logDenXParms<0)){
		logr = (scale*logDenPropParms - logDenXParms) / metropolisStepTemp["parms"]
		if ( is.numeric(logr) & (logr) <= log(runif(1)) ){
			#reject
			return(c(obs=NA, parms=-Inf))
		}
	}
	# evaluate the model at parameters theta0 
	tmp.pred <- fModel(theta0, ...)
	tmp.diffObs <- tmp.pred - obs
	logDenObs <- if( is.matrix(invCovar) ){
		t(tmp.diffObs) %*% invCovar %*% tmp.diffObs
	}else{
		sum(tmp.diffObs^2 / invCovar) # assuming invCovar specifies independent variances
	}
	tmp.misfit <-  c(obs=as.numeric(logDenObs), parms=as.numeric(logDenPropParms) )
	scale * tmp.misfit
	### the misfit: scale *( t(tmp.diffObs) %*% invCovar %*% tmp.diffObs + t(tmp.diffParms) %*% invCovarTheta %*% tmp.diffParms )
}
#mtrace(logDenGaussian)
#mtrace.off()
#twUtestF(logDenGaussian)

dummyTwDEMCModel <- function(
	### example model function: y=a+bx
	theta,	##<< parameter vector with names a and b
	xval	##<< additional argument, passed by ... in logDenGaussian
){ 
	# dummyTwDEMCModel
	##seealso<<   
	## \code{\link{twDEMCInt}}
	## \code{\link{logDenGaussian}}
	theta["a"] + theta["b"]*xval 
}
