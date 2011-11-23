twCalcLogDenPar <- function(
	### Invokes fLogDen with proposal in a parallel load balanced way.
	fLogDen,				##<< the objective function
	xProp,					##<< numeric matrix (nCases x nParm) of proposals 
	logDenCompX=NULL	
	### all components of logDensity of xProp (result of fLogDen * fLogDenScale) 
	### colnames must contain intResCompNames 
	### rows: number of cases in xProp	
	,intResCompNames=character(0)	
	### character vector: names of results components of fLogDen that are used for internal Metropolis decisions 
	,argsFLogDen=list()		##<< arguments passed to fLogDen
	,fLogDenScale=1			##<< factor multiplied to the result of fLogDen
	,debugSequential=FALSE	##<< see \code{\link{sfFArgsApplyLB}}
	,remoteDumpfileBasename=NULL,	##<< see \code{\link{sfRemoteWrapper}}
	...						##<< further arguments passed to fLogDen
){
	##seealso<< 
	## \code{\link{twDEMCInt}}
	#if( (0 == length(logDenCompX)) ) logDenCompX=xProp[,FALSE,drop=FALSE]
	if( {tmp<-list(...); any(""==names(tmp)) || length(names(tmp))!=length(tmp)} )
		("twCalcLogDenPar: encountered unnamed argument in ... Check for <- and ,, in list()")
	boProvideX2Argument <- (0 < length(intResCompNames))
	Lp <- fLogDenScale * if(boProvideX2Argument){
			#call fLogDen with second argument: the internal components
			if( 0 == length(logDenCompX) )
				logDenCompX <- matrix(-Inf, ncol=length(intResCompNames), nrow=nrow(xProp), dimnames=list(NULL,parms=intResCompNames))
			if( !is.numeric(logDenCompX) || !is.matrix(logDenCompX) || nrow(logDenCompX)!=nrow(xProp) )
				stop("logDenCompX must be a numeric matrix with one row for each chain and column names correponding to a subst of names of result vector of fLogDen")
			iNames <- match( intResCompNames, colnames(logDenCompX) )
			if( any(is.na(iNames)) )
				stop("if logDenCompX is given, it must contain named columns for each entry of intResCompNames")
			logDenCompXIntUnscaled <- logDenCompX[,iNames,drop=FALSE]/fLogDenScale	# maybe internally uses cost function -1/2*logDen instead of logDen	
			F_ARGS <- function(i){c(list(xProp[i,]),list(logDenCompXIntUnscaled[i,]))}
			#F_ARGS(1)
			resl <- sfFArgsApplyLB( nrow(xProp), F_ARGS, F_APPLY=sfRemoteWrapper, remoteFun=fLogDen	, debugSequential=debugSequential, remoteDumpfileBasename=remoteDumpfileBasename, SFFARGSAPPLY_ADDARGS=argsFLogDen, ...) 
			sfSimplifyLBResult(resl)
		}else{
			do.call( sfApplyMatrixLB, c(list( X=xProp, MARGIN=1, FUN=sfRemoteWrapper, remoteFun=fLogDen		, debugSequential=debugSequential, remoteDumpfileBasename=remoteDumpfileBasename), argsFLogDen, list(...)) )	#use doCall in order to use argsFLogDen
		}
	.logDen <- if( is.matrix(Lp) )
			colSums(Lp)
		else
			Lp
	.logDenComp <- if( is.matrix(Lp) )	t(Lp)	else matrix(Lp,ncol=1,dimnames=list(NULL,rownames(Lp)))
	list( logDen=.logDen, logDenComp=.logDenComp)
	### List with the following items \describe{
	### \item{logDen}{numeric vector: for each state: the sum of logDens over all components, multiplied by fLogDenScale}
	### \item{logDenComp}{numeric matrix: return components of fLogDen, one row for each state, columns: components }
	### }
}
attr(twCalcLogDenPar,"ex") <- function(){
	data(twdemcEx1)
	xProp <- stackChains(twdemcEx1$parms)
	
	data(twLinreg1)
	attach( twLinreg1 )
	res <- twCalcLogDenPar(logDenGaussian,xProp
		,fModel=dummyTwDEMCModel		### the model function, which predicts the output based on theta 
		,obs=obs			### vector of data to compare with
		,invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	str(res)
	resM <- matrix(res$logDen, ncol=dim(twdemcEx1$parms)[3])
	str(resM)
}

