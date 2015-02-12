
modTwTwoDenEx1 <- function(
        ### example model giving two predictions that can be compared to different observations
        theta	##<< model parameters a and b
        , xSparse	##<< numeric vector of Sparse input/output relationship 
        , xRich		##<< numeric vector of rich input/output relationship
        , thresholdCovar=0	##<< model structural deficiency
){
    ##details<< model output y1 represents a longterm observations
    ## It is based on longterm average of xRich instead of detailed values
    ## , Model output y1 represents a short measurement campaing. 
    ## During this campaing xSparse does not vary but detailed measurements of xRich are utilized 
    ## , In the short-term relation, the model may simulate a detailed threshold in the covariate
    ## or abstract from those details by thresholdCovar=0.
    
    ##value<< list with model predictions
    list( ##describe<<
            # with the following line, bias of theta[2] fully leaks into estimate of a
            #y1 = as.numeric(theta[1]*xSparse + theta[2]*mean(xRich))	##<< theta[1]*x1 + theta[2]*8
            y1 = as.numeric(theta[1]*xSparse + theta[2]*mean(xRich)/10)	##<< theta[1]*x1 + theta[2]*mean(xRich)/10
            ,y2 = as.numeric(theta[1]*xSparse[1] + theta[2]*pmax(0,xRich-thresholdCovar) ) 		 
    ) ##end<<
}




denSparse <- function(
	### Example of using two different logDensity functions: density of Sparse observations
	theta		
	,twTwoDenEx#=twTwoDenEx1
	,theta0=twTwoDenEx$thetaTrue
    ,intermediate=NULL
	,...
){
	theta0[names(theta)] <- theta
	#pred <- .memoize1( modTwoDenExCache, theta, { twTwoDenEx$fModel(theta0, xSparse=twTwoDenEx$xSparse, xRich=twTwoDenEx$xRich, ...) })
    pred <- twTwoDenEx$fModel(theta0, xSparse=twTwoDenEx$xSparse, xRich=twTwoDenEx$xRich, ...) 
    misfit <- twTwoDenEx$obs$y1 - pred$y1
	-1/2 * sum((misfit/twTwoDenEx$sdObs$y1)^2)
}


denSparsePrior <- function(
	### Same as denSpace but returning two components: prior and misfit
	theta		
	,twTwoDenEx#=twTwoDenEx1
	,theta0=twTwoDenEx$thetaTrue
	,thetaPrior = NULL		##<< the prior estimate of the parameters
	,invCovarTheta = NULL	##<< the inverse of the Covariance of the prior parameter estimates
    ,intermediate = list()  ##<< specify a mutable environement to store intermediate results that are the same across different densities for given theta
	,...
){
    ##details<< 
    ## Used for demonstrating usage of intermediate results.
    ## See test case \code{ofMultiIntermediate} in file unitTests/runittwDEMC.R
    #
    ##seealso<<
    ## \code{\link{denRichPrior}}
	theta0[names(theta)] <- theta
    logDenPropParms <- if( !is.null(thetaPrior) ){
                tmp.diffParms <- theta - thetaPrior
                if( is.matrix(invCovarTheta) )
                    t(tmp.diffParms) %*% invCovarTheta %*% tmp.diffParms
                else
                    sum(tmp.diffParms^2 / invCovarTheta )	# assume invCovar to be independent variances
            } else 0
    if( !length(intermediate) ){
        intermediate$pred <- twTwoDenEx$fModel(theta0, xSparse=twTwoDenEx$xSparse, xRich=twTwoDenEx$xRich, ...)
    }
    misfit <- twTwoDenEx$obs$y1 - intermediate$pred$y1
    ret <- -1/2 * c( parmsSparse=logDenPropParms , obsSparse=  sum((misfit/twTwoDenEx$sdObs$y1)^2) )
    attr(ret,"intermediate") <- intermediate
    ret
}


denRich <- function(
	### Example of using two different logDensity functions: density of data-rich observations
	theta
	,twTwoDenEx#=twTwoDenEx1
	,theta0=twTwoDenEx$thetaTrue
    ,modTwoDenExCache=NULL
    ,...
){
	theta0[names(theta)] <- theta
	#pred <- .memoize1( modTwoDenExCache, theta, {twTwoDenEx$fModel(theta0, xSparse=twTwoDenEx$xSparse, xRich=twTwoDenEx$xRich,...)})
    pred <- twTwoDenEx$fModel(theta0, xSparse=twTwoDenEx$xSparse, xRich=twTwoDenEx$xRich,...)
    misfit <- twTwoDenEx$obs$y2 - pred$y2
	-1/2 * sum((misfit/twTwoDenEx$sdObs$y2)^2)
}

denRichPrior <- function(
        ### Example of using two different logDensity functions: density of data-rich observations
        theta
        ,twTwoDenEx#=twTwoDenEx1
        ,theta0=twTwoDenEx$thetaTrue
        ,thetaPrior = NULL		##<< the prior estimate of the parameters
        ,invCovarTheta = NULL	##<< the inverse of the Covariance of the prior parameter estimates, or alternatively the diag of Covariance matrix (sigma_i)
        ,intermediate = list()  ##<< specify a mutable environement to store intermediate results that are the same across different densities for given theta
        ,...
){
    ##seealso<<
    ## \code{\link{denSparsePrior}}
    #
    theta0[names(theta)] <- theta
    logDenPropParms <- if( !is.null(thetaPrior) ){
                tmp.diffParms <- theta - thetaPrior
                if( is.matrix(invCovarTheta) )
                    t(tmp.diffParms) %*% invCovarTheta %*% tmp.diffParms
                else
                    sum(tmp.diffParms^2 / invCovarTheta )	# assume invCovar to be independent variances
            } else 0
    if( !length(intermediate) ){
        intermediate$pred <- twTwoDenEx$fModel(theta0, xSparse=twTwoDenEx$xSparse, xRich=twTwoDenEx$xRich, ...)
    } 
    misfit <- twTwoDenEx$obs$y2 - intermediate$pred$y2
    ret <- -1/2 * c( parmsRich=logDenPropParms , obsRich=  sum((misfit/twTwoDenEx$sdObs$y2)^2) )
    attr(ret,"intermediate") <- intermediate
    ret
}

denBoth <- function(
	### Example of using two different logDensity functions: comparison of combining denSparse and denRich into one function.
	theta
	,twTwoDenEx#=twTwoDenEx1
	,weights=c(1,1,1)			# weights for the two data streams
	,theta0=twTwoDenEx$thetaTrue
	,thetaPrior = NULL		##<< the prior estimate of the parameters
	,invCovarTheta = NULL	##<< the inverse of the Covariance of the prior parameter estimates
    ,modTwoDenExCache=NULL          ##<< environment to cache results of model prediction for given theta
    ,sdSparse=twTwoDenEx$sdObs$y1   ##<< standard deviation of sparse observations
    ,sdRich=twTwoDenEx$sdObs$y2     ##<< standard deviation of rich observations
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
	#pred <- .memoize1( modTwoDenExCache, theta, {twTwoDenEx$fModel(theta0, xSparse=twTwoDenEx$xSparse, xRich=twTwoDenEx$xRich, ...)})
    pred <- twTwoDenEx$fModel(theta0, xSparse=twTwoDenEx$xSparse, xRich=twTwoDenEx$xRich, ...)
    misfit <- twTwoDenEx$obs$y1 - pred$y1
	d1 <- sum((misfit/sdSparse)^2)
	misfit <- twTwoDenEx$obs$y2 - pred$y2
	d2 <- sum((misfit/sdRich)^2)
	db <- -1/2 * c(parms=logDenPropParms, y1=d1,y2=d2)*weights
	db
}
attr(denBoth,"ex") <- function(){
	data(twTwoDenEx1)
	denBoth(twTwoDenEx1$thetaTrue, twTwoDenEx1)
	denBoth(twTwoDenEx1$thetaTrue, twTwoDenEx1, thresholdCovar = 0.3)
}
