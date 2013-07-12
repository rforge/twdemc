twCalcLogDenPar <- function(
	### Invokes fLogDen with proposal in a parallel load balanced way.
	fLogDen,				##<< the objective function
	xProp,					##<< numeric matrix (nCases x nParm) of proposals 
	logDenCompX=NULL	
	### all components of logDensity of xProp (result of fLogDen) 
	### colnames must contain intResCompNames 
	### rows: number of cases in xProp	
	,intResComp=character(0)	
	### character vector: names of results components of fLogDen that are used for internal Metropolis decisions 
	,argsFLogDen=list()		##<< arguments passed to fLogDen
	,debugSequential=FALSE	##<< see \code{\link{sfFArgsApplyLB}}
	,remoteDumpfileBasename=NULL,	##<< see \code{\link{sfRemoteWrapper}}
	...						##<< further arguments passed to fLogDen
){
	##seealso<< 
	## \code{\link{twDEMCBlockInt}}
	#if( (0 == length(logDenCompX)) ) logDenCompX=xProp[,FALSE,drop=FALSE]
	if( {tmp<-list(...); any(""==names(tmp)) || length(names(tmp))!=length(tmp)} )
		("twCalcLogDenPar: encountered unnamed argument in ... Check for <- and ,, in list()")
	boProvideX2Argument <- (0 < length(intResComp))
	Lp <- if(boProvideX2Argument){
			#call fLogDen with second argument: the internal components
			if( 0 == length(logDenCompX) )
				logDenCompX <- matrix(-Inf, ncol=length(intResComp), nrow=nrow(xProp), dimnames=list(NULL,parms=intResComp))
			if( !is.numeric(logDenCompX) || !is.matrix(logDenCompX) || nrow(logDenCompX)!=nrow(xProp) )
				stop("logDenCompX must be a numeric matrix with one row for each chain and column names correponding to a subst of names of result vector of fLogDen")
			iNames <- match( intResComp, colnames(logDenCompX) )
			if( any(is.na(iNames)) )
				stop("if logDenCompX is given, it must contain named columns for each entry of intResCompNames")
			logDenCompXIntUnscaled <- logDenCompX[,iNames,drop=FALSE]		
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
	### \item{logDen}{numeric vector: for each state: the sum of logDens over all components}
	### \item{logDenComp}{numeric matrix: return components of fLogDen, one row for each state, columns: components }
	### }
}
attr(twCalcLogDenPar,"ex") <- function(){
	data(twdemcEx1)
	xProp <- stackChains(twdemcEx1)[,-1]
	
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
}

.twCalcLogDens <- function( 
     ### calculating all densities for given proposal (called remotely from \code{\link{twCalcLogDensPar}})
     x				    ##<< numeric vector (nParm)
    ,dInfos				##<< list describing the logDensities, see \code{\link{twDEMCBlockInt}}
){
    intermediates <- list()
    #recover()
    #dInfo <- dInfos[[1]]
    #dInfo <- dInfos[[2]]
    #dInfo <- dInfos[[3]]
    LpL <- list()
    for( dInfo in dInfos){
        iName <- dInfo$intermediate
        if( !is.null(iName) ) dInfo$argsFLogDen$intermediate <- intermediates[[iName]]       # set the intermediate
        Lp <- do.call( dInfo$fLogDen, c( list(x), dInfo$argsFLogDen ))
        if( !is.null(iName) ) intermediates[[ iName ]] <- attr(Lp, "intermediate")
        LpL[[ length(LpL)+1 ]] <- structure( as.numeric(Lp), names=names(Lp) )  # drop all other attributes on storing
    }
    #names(LpL) <- NULL  # prevent recreating names in c
    # need to return the full list instead of vector, because number of components is important
    ### List with entry for ech logDensity: numeric vector: result components of the logDensity
    LpL    
}


twCalcLogDensPar <- function(
        ### Invokes all fLogDens with proposal in a parallel way, taking care of intermediates between densities
        dInfos				##<< list describing the logDensities, see \code{\link{twDEMCBlockInt}}
        ,xProp				##<< numeric matrix (nCases x nParm) of proposals 
        ,debugSequential=FALSE	##<< see \code{\link{sfFArgsApplyLB}}
        ,remoteDumpfileBasename=NULL	##<< see \code{\link{sfRemoteWrapper}}
){
    ##seealso<< 
    ## \code{\link{twDEMCBlockInt}}, \code{\link{twCalcLogDenPar}}
    #
    ##details<<
    ## Does not take care of internal components: provides no argument logDenPrev
    #
    if( nrow(xProp)==0 ){
        stop("twCalcLogDensPar: encountered empty parameter matrix. Cannot infor result component names.")
    } else if( nrow(xProp)==1 ){
        # degenrate case of only one row
        x <- xProp[1,]
        LpVecL <- .twCalcLogDens( x, dInfos=dInfos )
        LpVec <- do.call(c,LpVecL)
        Lp <- matrix( LpVec, nrow=1, dimnames=list(NULL, names(LpVec) ) )
        .logDen <- sum(LpVec)
    }else{
        # more than one row
        Lpl <- if( debugSequential || !sfParallel() ){
            Lpl <- (apply( xProp, 1, .twCalcLogDens, dInfos=dInfos ))
            #row <- Lpl[[1]]
        }else{
            if( length(remoteDumpfileBasename) ){
                # for debugging use RemoteWrapper to write dumpfile
                Lpl <- (sfApply( xProp,1, sfRemoteWrapper, remoteFun=.twCalcLogDens
                                , dInfos = dInfos
                                , remoteDumpfileBasename=remoteDumpfileBasename )
                )
            }else{
                # no need of remote wrapper, which is faster 
                Lpl <- (sfApply( xProp,1, .twCalcLogDens, dInfos = dInfos ))
            }
        }
        LpVecL <- Lpl[[1]]
        Lpt <- sapply( Lpl, function(row){ do.call(c,row)} )
        if( !is.matrix(Lpt) ){
            # degenerate case of returning only one component, hence Lpt is a vector
            Lp <-  matrix(Lpt, ncol=1, dimnames=list(NULL,names(dInfos)) )
            .logDen <- Lpt
        } else{
            # several rows
            Lp <- t(Lpt) 
            .logDen <- rowSums(Lp)
        } 
    }
    #LpVecL is one result of call to .twCalcLogDens a list with components  
    LpPos <- rep(seq_along(LpVecL), sapply(LpVecL,length) )
    list( logDen=.logDen, logDenComp=Lp, logDenCompPos=LpPos)
    ### List with the following items \describe{
    ### \item{logDen}{numeric vector: for each state: the sum of logDens over all components}
    ### \item{logDenComp}{numeric matrix: return components of fLogDen, one row for each state, columns: components }
    ### \item{logDenCompPos}{integer vector (nDenComp): index of the densitiy that provides the component }
    ### }
}
attr(twCalcLogDensPar,"ex") <- function(){
    # see test function twCalcLogDensPar in test case twCalcLogDenPar
}


