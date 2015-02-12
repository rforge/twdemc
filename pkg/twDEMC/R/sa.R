
twDEMCSA <- function(
	### Simulated annealing DEMC 
	thetaPrior			##<< vector of parameters, point estimate
        ##, alternatively array with initial states, as returned by \code{\link{initZtwDEMCNormal}}
	,covarTheta			##<< the a prior covariance of parameters, see \code{link{initZtwDEMCNormal}}
	,nGen=512			##<< number of generations, i.e. Monte Carlo steps
	,nObs				##<< integer vector specifying the number of observations for each result component
        ## for each name in result of logDensity functions there must be an entry in nObs
	,... 				##<< further argument to \code{\link{twDEMCBlockInt}}
	, dInfos			##<< Information on density functions (see \code{\link{twDEMCBlockInt}})
	, m0 = calcM0twDEMC(length(thetaPrior),nChainPop)	##<< minimum number of samples in step for extending runs
	, controlTwDEMC = list()	##<< control parameters with entry \code{thin} (see \code{\link{twDEMCBlockInt}})
	, debugSequential=FALSE		##<< set to TRUE to avoid parallel execution, good for debugging
	, restartFilename=NULL		##<< filename to write intermediate results to
    , remoteDumpfileBasename=NULL		##<< fileBasename to write dumps to on error
    #
	,nChainPop=4		##<< number of chains within population
	,nPop=2				##<< number of populations
	,doIncludePrior=FALSE	##<< should the prior be part of initial population 
		##<< Recommendation to set to false, because if the TRUE parameter is in initial set, initial Temperature calculation is too optimistic
	#
	, ctrlBatch = list(             ##<< list of arguments controlling batch executions, see also \code{\link{twDEMCSACont}} 
		##describe<< 
		nGen0=m0*controlTwDEMC$thin*3	##<< number of generations for the initial batch
		#,useSubspaceAdaptation=FALSE	# #<< if TRUE then overall space is devided and each subspace is explored with locally adapted DEMC, see	\code{\link{divideTwDEMCSACont}}	
		##end<<
		)	
	, ctrlT = list(    ##<< list of arguments controlling Temperature decrease, see also \code{\link{twDEMCSACont}} 
		##describe<< 
		qTempInit=0.4		##<< quantile of logDensities used to calculate initial beginning and end temperature, with default 0.4: 40% of the space is accepted
        ,TBaseInit=NULL     ##<< numeric scalar: initial base temperature. If given this is used for calculating initial temperature
        ,isVerbose=FALSE    ##<< boolean scalar: set to TRUE to report stream base temperatures during batches
        ,TEndFixed=NULL     ##<< set to a scalar end temperature, e.g. in order to decrease temperatue to a given temperature
        ,qBest=0.1          ##<< the proportion of best samples to base calculation of target temperature on 
        ##end<<
	)
	, ctrlConvergence = list()		##<< list or arguments controlling check for convergence, see \code{\link{twDEMCSACont}}  
	#, ctrlSubspaces = list()		# #<< list of arguments controlling splitting and merging of subspaces, see \code{\link{divideTwDEMCSACont}}
){
	##detail<< \describe{\item{Continuing a previous run}{
	## When supplying the first argument an object of class twDEMCPops to twDEMCSA, then this run is extended.
	## All other parameters, except nGen, are then are ignored.
	## In order to change parameters, modify list entry args in the twDEMCPops object.
	##}}
	if( inherits(thetaPrior,"twDEMCPops") && 0 != length(thetaPrior$args) ){
		if( !missing(debugSequential) )
			thetaPrior$args$debugSequential=debugSequential
        if( !missing(restartFilename) )
            thetaPrior$args$restartFilename=restartFilename
        if( !missing(remoteDumpfileBasename) )
            thetaPrior$args$remoteDumpfileBasename=remoteDumpfileBasename
        if( !missing(controlTwDEMC) )
            thetaPrior$args$controlTwDEMC <- twMergeLists( thetaPrior$args$controlTwDEMC, controlTwDEMC )
        if( !missing(ctrlBatch) )
            thetaPrior$args$ctrlBatch <- twMergeLists( thetaPrior$args$ctrlBatch, ctrlBatch )
        if( !missing(ctrlT) )
            thetaPrior$args$ctrlT <- twMergeLists( thetaPrior$args$ctrlT, ctrlT )
        if( !missing(ctrlConvergence) )
            thetaPrior$args$ctrlConvergence <- twMergeLists( thetaPrior$args$ctrlConvergence, ctrlConvergence )
        #if( !missing(ctrlSubspaces) )
        #    thetaPrior$args$ctrlSubspaces <- twMergeLists( thetaPrior$args$ctrlSubspaces, ctrlSubspaces )
        if( !missing(dInfos) ){
            warning("twDEMCSA: continueing previous run, but with different argument dInfos")
            thetaPrior$args$dInfos <- dInfos   # overwrite
        }
        thetaPrior$failureMsg <- NULL       # clear previous return and error messages
        thetaPrior$finishEarly <- NULL
        #ret <- if( 0!=length(thetaPrior$args$ctrlBatch) &&
		#		   0!=length(thetaPrior$args$ctrlBatch$useSubspaceAdaptation) &&	
		#		   thetaPrior$args$ctrlBatch$useSubspaceAdaptation 
		#	)
		#	ret <- do.call(divideTwDEMCSACont,c(list(thetaPrior, nGen=nGen),thetaPrior$args))
		#else 	
		mcPops <- do.call(twDEMCSACont,c(list(thetaPrior, nGen=nGen),thetaPrior$args))
		#ret$args$ctrlBatch$useSubspaceAdaptation <- ctrlBatch$useSubspaceAdaptation
		return(mcPops)
	}
    #
    #-------------- generate initial sample
    # check for array before evaluating formal arguments, because m0 then changes, which is used in default values
    #mtrace(initZtwDEMCNormal)
    Zinit0 <- if( is.array(thetaPrior) ){
                .dimZ <- dim(thetaPrior)
                if( length(.dimZ) != 3 ) stop("twDEMCSA: wrong dimensions of thetaPrior")
                nChain <- .dimZ[3] 
                if( missing(nChainPop) ) nChainPop <- nChain %/% nPop
                if( missing(nPop) ) nPop <- nChain %/% nChainPop
                if( .dimZ[3] != nChainPop*nPop ) stop("twDEMCSA: wrong number of chains in thetaPrior")
                m0 = calcM0twDEMC( ncol(thetaPrior) ,nChainPop)	##<< minimum number of samples in step for extending runs
                if( .dimZ[1] < m0 ) stop(paste("twDEMCSA: too few cases in  thetaPrior (",.dimZ[1],"). Need at least m0=",m0,sep=""))
                thetaPrior
            } else initZtwDEMCNormal( thetaPrior, covarTheta, nChainPop=nChainPop, nPop=nPop, doIncludePrior=doIncludePrior, m0=m0)
    #
	#--------  fill in default argument values
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
	thin <- controlTwDEMC$thin
    # The user may provide list arguments to overwrite some default entries.
    # With formals and twMergeLists, we make sure to take the other default entries in the lists
	frm <- formals()
	#ctrlSubspaces <- if( hasArg(ctrlSubspaces) ) twMergeLists( eval(frm[["ctrlConvergence"]]), ctrlSubspaces ) else ctrlSubspaces
	ctrlBatch <- if( hasArg(ctrlBatch) ) twMergeLists( eval(frm[["ctrlBatch"]]), ctrlBatch ) else ctrlBatch	
	ctrlT <- if( hasArg(ctrlT) ) twMergeLists( eval(frm[["ctrlT"]]), ctrlT ) else ctrlT
    #
    #---------- calculate logDensity of initial components
	ss <- stackChains(Zinit0)
    #print("before twCalcLogDensPar"); recover()
    cat("Calculating logDensity for ", length(ss)," initial states for m0=",m0,"rows per population.\n")
    logDenDS0 <- (tmp <- twCalcLogDensPar(dInfos, ss, remoteDumpfileBasename=remoteDumpfileBasename))$logDenComp
	boFinite <- apply(is.finite(logDenDS0), 1, all )
	m0FiniteFac <- sum(boFinite) / nrow(logDenDS0)
    # if too few initial states have yielded finite results of logDensity functions, add some more initial states
	if( m0FiniteFac == 1){
			Zinit <- Zinit0
			logDenDS <- logDenDS0
		}else if( m0FiniteFac < 0.05 ){
			stop(paste("twDEMCSA: less than 5% finite solutions in initial exploration (",m0FiniteFac,"). Check setup."),sep="")
		}else if( m0FiniteFac > 0.9 ){
			sumLogDen0 <- rowSums(logDenDS0)
			Zinit <- replaceZinitNonFiniteLogDens( Zinit0, sumLogDen0 )
			logDenDS <- logDenDS0
			logDenDS[!boFinite,] <- logDenDS[which.min(sumLogDen0),]
		}else{
			# generate more proposals and concatenate finite cases from both
			Zinit1 <- initZtwDEMCNormal( thetaPrior, covarTheta, nChainPop=nChainPop, nPop=nPop, doIncludePrior=doIncludePrior
				,m0FiniteFac=min(1,m0FiniteFac)		
			)
			ss1 <- stackChains(Zinit1)
			logDenL <- lapply( dInfos, function(dInfo){
					resLogDen <- twCalcLogDenPar( dInfo$fLogDen, ss1, argsFLogDen=dInfo$argsFLogDen)$logDenComp
				})
			# replace missing cases
			logDenDS1 <- abind( logDenL )
			# logDenDS1[ 1:round(nrow(logDenDS1)/1.7),1] <- NA		# for testing replacement of non-finite cases
			boFinite1 <- apply(is.finite(logDenDS1), 1, all )
			m0FiniteFac1 <- sum(boFinite1) / nrow(logDenDS1)
			ss12 <- rbind( ss[boFinite, ,drop=FALSE], ss1[boFinite1, ,drop=FALSE] )
			logDenDS12 <- rbind( logDenDS0[boFinite, ,drop=FALSE], logDenDS1[boFinite1, ,drop=FALSE])
			nDiff <- nrow(ss) - nrow(ss12)
			if( nDiff/nrow(ss) > 0.1 ){
				stop("twDEMCSA: too many states yielding non-finite logDensities.")
			}else if( nDiff > 0 ){
				# duplicate missings to refill
				iSample <- sample( nrow(ss12), size=nDiff )
				ss12 <- rbind( ss12, ss12[ iSample, ] )
				logDenDS12 <- rbind( logDenDS12, logDenDS12[iSample,] )
			}else if(nDiff < 0){
				iSample <- sample( nrow(ss12), size=nrow(ss))
				ss12 <- ss12[ iSample, ]
				logDenDS12 <- logDenDS12[iSample,]
			}  
			## reshape to chains
			ncolZ <- dim(Zinit0)[3]
			tmp <- matrix(1:nrow(ss12), ncol=ncolZ)
			Zinit <- abind( lapply( 1:ncolZ, function(i){ ss12[ tmp[,i], ,drop=FALSE]}), rev.along=0 )
			logDenDS <- logDenDS12
	} # end generating nonfinite Zinit
    #
    # apply overfitting correction?
    # No need here because Temperature would be set to 1 for corrected and noncorrected values
    #       
	# now the length of resComp is known (colnames(logDenDS)), so check the TFix and TMax parameters
    iNoNames <- which( "" == colnames(logDenDS) )
    if( length(iNoNames) ) stop("Encountered unnamed logDensitiy components. Provide names of each return component in fLogDen.") 
    nResComp <- ncol(logDenDS)
	ctrlT$TFix <- completeResCompVec( ctrlT$TFix, colnames(logDenDS) )
	iFixTemp <- which( is.finite(ctrlT$TFix) )
    iNonFixTemp <- which( !is.finite(ctrlT$TFix) )
	ctrlT$TMax <- tmp <- completeResCompVec( ctrlT$TMax, colnames(logDenDS) )
	iMaxTemp <- which( is.finite(ctrlT$TMax) )	
    #
    # check nObs names and set all TFix to 1 
    if( 0 == length(names(nObs))){
        if( length(nObs) == length(iNonFixTemp) ) names(nObs) <- colnames(logDenDS)[iNonFixTemp] else
           stop("Provide names with nObs, or provide one component for each temperated (i.e. not in TFix) fLogDen return component.")
    }
    iMissing <- which( !(colnames(logDenDS)[iNonFixTemp] %in% names(nObs)) )
    if( length(iMissing) ) stop("Missing the following components in nObs: ",paste( colnames(logDenDS)[iNonFixTemp][iMissing],sep=",") )
    # make nObs correspond to each data stream, set TFix Components to 1    
    nObsComp <- structure( rep( NA_real_, ncol(logDenDS) ), names=colnames(logDenDS) )
    #nObsComp[names(nObs)] <- nObs  # if there are seveval occurence of the same name, others still NA
    # name <- names(nObs)[1]
    for( name in names(nObs) ){  nObsComp[ names(nObsComp) == name ] <- nObs[name]    }
    nObs <- nObsComp
	#
	#if( 0 == length( TMaxInc) ) TMaxInc <- structure( rep( NA_real_, nResComp), names=colnames(logDenDS) )
	#if( nResComp != length( TMaxInc) ) stop("twDEMCSA: TMaxInc must be of the same length as number of result Components.")
	#iMaxIncTemp <- which( is.finite(TMaxInc) )
	#
	#------ intial temperatures: 
	#sapply( seq_along(expLogDenBest), function(i){ max(logDenDS[,i]) })
	##details<< \describe{\item{Initial temperature}{
	## Initial parameters are ranked according to their maximum log-Density across components.
	## The parameters and logDensity results at rank position defined by argument \code{ctrlT$qTempInit} is selected.
	## The stream temperatures are inferred by deviding logDensity components by the number of observations.
    ## From these, a common base temperatue is calculated by (see \code{\link{calcBaseTemp}}) and rescaled to stream temperatures by \code{\link{calcStreamTemp}}.
	##}}
    #print("twDEMCSA: before calculating initial temperature"); recover()
    if( length(ctrlT$TBaseInit) ){
        temp0 <- calcStreamTemp(ctrlT$TBaseInit, nObs, TFix=ctrlT$TFix, iFixTemp=iFixTemp)        
    } else {
        iTmax <- 5
        T0s <- rep(NA_real_, iTmax )
        T <- nObs
        T[iFixTemp] <- ctrlT$TFix[iFixTemp] 
        T0s[1] <- .Machine$double.xmax
        relChange <- relChangeReq <- 0.05   # if below 5% change then no need to iterate further
        iT <- 1
        while( (iT < iTmax) && (abs(relChange) >= relChangeReq) ){
            iT <- iT + 1
            resLogDenT <- calcTemperatedLogDen( logDenDS, T )
            iBest <- getBestModelIndices( resLogDenT, dInfos, prob=ctrlT$qTempInit )
            T0 <- T0s[iT] <- calcBaseTempSk( -2*logDenDS[iBest, ,drop=FALSE], nObs, ctrlT$TFix, iFixTemp = iFixTemp, iNonFixTemp = iNonFixTemp)
            T <- calcStreamTemp(T0, nObs, ctrlT$TFix, iFixTemp = iFixTemp)
            relChange <- (T0 - T0s[iT-1])/T0 
        }
        temp0 <- T
    	if( !all(is.finite(temp0)) ) stop("twDEMCSA: encountered non-finite Temperatures.")
    }
    ctrlT$TMax[iNonFixTemp] <- pmin(ifelse(is.finite(ctrlT$TMax),ctrlT$TMax, Inf), ifelse(is.finite(temp0),temp0,Inf) )[iNonFixTemp]		# decrease TMax  
    #if( any(temp0 > 8)) stop("twDEMCSA: encountered too high temperature.")	
    print(paste("initial T=",paste(signif(temp0,2),collapse=",")," for ",min(nGen, ctrlBatch$nGen0)," generations","    ", date(), sep="") )
	#if( sum(temp0 > 1) == 0){ dump.frames("parms/debugDump",TRUE); stop("twDEMCSA: no temperature > 0") }	
	#
	mcPops <- res0 <-  res <- twDEMCBlock( Zinit
		, nGen=min(nGen, ctrlBatch$nGen0)
		#,nGen=16, debugSequential=TRUE
		, dInfos=dInfos
		, TSpec=cbind( T0=temp0, TEnd=temp0 )
		, nPop=nPop
		, m0=m0, controlTwDEMC=controlTwDEMC, debugSequential=debugSequential
		# XXTODO: replace lines below later on by ...
		#, blocks = blocks
		,...
	)
    print(paste("finished initial ",ctrlBatch$nGen0," out of ",nGen," gens.    ", date(), sep="") )
	#
    #if( nGen > ctrlBatch$nGen0 ){ #depre, now need to call Continue to return proper diagnostics
    #if( TRUE ){     # need to call Continue to return proper diagnostics
		#nObsLocal <- nObs 
		#ret <- if( ctrlBatch$useSubspaceAdaptation ){
		#	ret <- divideTwDEMCSACont( mc=res0, nGen=nGen-ctrlBatch$nGen0, nObs=nObs
		#		, m0=m0, controlTwDEMC=controlTwDEMC, debugSequential=debugSequential, restartFilename=restartFilename
		#		, ctrlBatch=ctrlBatch, ctrlT=ctrlT, ctrlConvergence=ctrlConvergence, ctrlSubspaces=ctrlSubspaces
		#		, ...
		#	) 
		#}else {
			mcPops <- twDEMCSACont( mc=res0, nGen=nGen-ctrlBatch$nGen0, nObs=nObs
				, m0=m0, controlTwDEMC=controlTwDEMC, debugSequential=debugSequential, restartFilename=restartFilename
				, ctrlBatch=ctrlBatch, ctrlT=ctrlT, ctrlConvergence=ctrlConvergence
				, ...
				)
		#}
	#}	
	#ret$args$ctrlBatch$useSubspaceAdaptation <- ctrlBatch$useSubspaceAdaptation	# remember this argument
    ##value<< An object of class \code{twDEMCPops} as described in  \code{\link{twDEMCBlockInt}}.
    ## See \code{\link{subset.twDEMCPops}} for processing further handling of this class.
	mcPops
}
attr(twDEMCSA,"ex") <- function(){
# Best to have a look at the the package vignettes
    
#--------------- single density ----------------------------
# We will use logDenGaussian as logDensity function that compares a simple linear model with observations.
data(twLinreg1)

# collect all the arguments to the logDensity in a list (except the first argument of changing parameters)    
argsFLogDen <- with( twLinreg1, list(
        fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
        obs=obs,			    ### vector of data to compare with
        invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
        thetaPrior = thetaTrue,	### the prior estimate of the parameters
        invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
        xval=xval
))
do.call( logDenGaussian, c(list(theta=twLinreg1$theta0),argsFLogDen))
do.call( logDenGaussian, c(list(theta=twLinreg1$thetaTrue),argsFLogDen))    # slightly largere misfit than nObs/2=15, underestimated sdObs

.nGen=200
.nPop=2
mcPops <-  twDEMCSA( 
        theta=twLinreg1$theta0, covarTheta=diag(twLinreg1$sdTheta^2)       # to generate initial population
        , nGen=.nGen
        , dInfos=list(den1=list(fLogDen=logDenGaussian, argsFLogDen=argsFLogDen))
        , nPop=.nPop                                        # number of independent populations
        , controlTwDEMC=list(thin=4)                        # see twDEMCBlockInt for all the tuning options
        , ctrlConvergence=list(maxRelTChangeCrit=0.1)       # ok if T changes less than 10% 
        , ctrlT=list(TFix=c(parms=1))                       # do not use increased temperature for priors
        , nObs=c(obs=length(argsFLogDen$obs))               # number of observations used in temperature calculation
)
#mcp <- twDEMCSA( mcp, nGen=2000) 
mcPops <- twDEMCSA( mcPops, nGen=400)     # continue run 

rescoda <- as.mcmc.list(mcPops)
plot(rescoda, smooth=FALSE)
mcChains1 <- concatPops(mcPops)                   # array representation instead of list of pops, last dim is the chain
mcChains2 <- concatPops(stackChainsPop(mcPops))   # combining dependent chains within one population
mcChains3 <- concatPops(subsetTail(mcPops,0.5))   # take only the last part of the chains
c(getNGen(mcChains1), getNGen(mcChains3))
plot(as.mcmc.list(mcChains2), smooth=FALSE)



#--------------- multiple densities -------------------------
data(twTwoDenEx1)

thetaPrior <- twTwoDenEx1$thetaTrue
covarTheta <- diag((thetaPrior*0.3)^2)
invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation

thresholdCovar = 0.3	# the true value used to generate the observations
thresholdCovar = 0		# the effective model that glosses over this threshold

#str(twTwoDenEx1)
nObs <- c( obsSparse=length(twTwoDenEx1$obs$y1), obsRich=length(twTwoDenEx1$obs$y2) )

dInfos=list(
	dSparse=list(fLogDen=denSparsePrior
        , argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta)
        #, maxLogDen=-1/2*nObs[c("parmsSparse","obsSparse")] # control overfitting
        )
	,dRich=list(fLogDen=denRichPrior
        , argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta)
        , maxLogDen=c(parmsRich=0,-1/2*nObs["obsRich"])     # control overfitting for rich datastream
        )
)
blocks = list(
	a=list(dInfoPos="dSparse", compPos="a")
	,b=list(dInfoPos="dRich", compPos="b")
)

names(do.call( dInfos$dSparse$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparse$argsFLogDen)))
names(do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen)))


#trace(twDEMCSACont, recover )
#trace(twDEMCSA, recover )
res <- res0 <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs
	, nGen=3*256
	, ctrlT=list( TFix=c(parmsSparse=1,parmsRich=1) )   # no increased Temperature for priors
	, ctrlBatch=list( nGenBatch=256 )
	, debugSequential=TRUE
    , controlTwDEMC = list(
           DRgamma=0.1                          # use Delayed rejection
           #,controlOverfittingMinNObs = 20      # use overfitting control (for obsRich), recommended on using single density 
    )
	#, restartFilename=file.path("tmp","example_twDEMCSA.RData")
)
res <- twDEMCSA( res0, nGen=2*256 )	# extend the former run

(TCurr <- getCurrentTemp(res))
# Note that T does decrease to 1
# This accounts for structural model mismatch in addition to observation uncertinaty

mc0 <- concatPops(res)
mcE <- concatPops(subsetTail(res,0.2))      # only the last 20%
plot( as.mcmc.list(mc0) , smooth=FALSE )
matplot( mc0$temp, type="l" )
logDenT <- calcTemperatedLogDen(stackChains(mcE$resLogDen), TCurr)
iBest <- getBestModelIndex( logDenT, res$dInfos )
maxLogDenT <- logDenT[iBest, ]
ss <- stackChains(mcE$parms)
(thetaBest <- ss[iBest, ])
twTwoDenEx1$thetaTrue
(.qq <- apply(ss,2,quantile, probs=c(0.025,0.5,0.975) ))    # model error really in paramter b
plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twMisc::twRescale(rowSums(logDenT),c(10,100))] )
plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twMisc::twRescale(logDenT[,"obsSparse"],c(10,100))] )
plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twMisc::twRescale(logDenT[,"obsRich"],c(10,100))] )
plot( ss[,"a"], ss[,"b"], col=rgb(
		twMisc::twRescale(logDenT[,"obsSparse"]),0, twMisc::twRescale(logDenT[,"obsRich"]) ))
apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )

# density of parameters
plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")

# predictive posterior (best model only)
pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparse=xSparse, xRich=xRich) )
plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1)    # note that deviation is now in y2 - consistent with the introduced bias
}       # end examples twDEMCSA

#plot predictive posterior
.tmp.f <- function(){
    ss <- stackChains(d3$parms)
    apply( ss, 2, mean )
    apply( ss, 2, sd )
    nSamp <- 40
    pSamp <- ss[ sample.int(nSamp), ]
    pPred <- apply( pSamp, 1, fModel, xval=xval )
    statPred <- t(apply( pPred, 1, function(obsi){ c(mean(obsi), sd(obsi) )} ))
    plot(statPred[,1] ~ xval, ylim=range(statPred[,1] +c(1.96,-1.96)*statPred[,2]))
    lines( statPred[,1] +1.96*statPred[,2] ~ xval )
    lines( statPred[,1] -1.96*statPred[,2] ~ xval )
    points( obs ~ xval, col="maroon")
    
    w <- 1/sdObs^2/sum(1/sdObs^2)
    lm2 <- lm( obs ~ xval, weights=w )
    w <- 1/sdObs^2      # would more correspond to real case see ?lm, but summed to 1 internally
    lm2 <- lm.wfit( cbind(1,xval), obs, w=w )
    summary(lm2)
    
    obsTrue <- fModel( thetaTrue, xval)
    plot(statPred[,1] ~ obsTrue, ylim=range(statPred[,1] +c(1.96,-1.96)*statPred[,2]))
    lines( statPred[,1] +1.96*statPred[,2] ~ obsTrue )
    lines( statPred[,1] -1.96*statPred[,2] ~ obsTrue )
    lines( statPred[,1] -1.96*sdObs ~ obsTrue, col="maroon" )
    lines( statPred[,1] +1.96*sdObs ~ obsTrue, col="maroon" )
    abline(0,1)
    points( obs ~ obsTrue, col="maroon")
    
}


twDEMCSACont <- function(
	### continuing simulated annealing DEMC based on previous result
	mc					##<< result of twDEMCBlock
	,nGen=512			##<< overall number of generations to add (within one batch provide ctrlBatch$nGenBatch)
	,nObs				##<< integer vector (nResComp) specifying the number of observations for each result component
	,... 				##<< further argument to \code{\link{twDEMCBlockInt}}
	, m0 = calcM0twDEMC(getNParms(mc),getNChainsPop(mc))	##<< minimum number of samples in step for extending runs
	, controlTwDEMC = list()		##<< list argument to \code{\link{twDEMCBlockInt}} containing entry thin
	, debugSequential=FALSE		##<< set to TRUE to avoid parallel execution, good for debugging
	, restartFilename=NULL		##<< filename to write intermediate results to
	#
	, ctrlBatch = list(				##<< list of arguments controlling batch executions
		##describe<< 
        nSamplesBatch=m0*8          ##<< number of samples for one call to twDEMCStep (multiplied by thin to get generations), defaults to 8*m0
		#nGenBatch=m0*controlTwDEMC$thin*8		##<< number of generations for one call to twDEMCStep
        #nGenBatch=m0*controlTwDEMC$thin*3		##<< number of generations for one call to twDEMCStep
        ##<< default: set in a way that on average each population (assuming half are significant) is appended by 2*m0 samples
		#, nSampleMin=32				##<< minimum number of samples in each population within batch so that calculation of average density is stable
		, pThinPast=0.5				##<< in each batch thin the past to given fraction before appending new results
		##end<<
	)	
	, ctrlT = list(					##<< list of arguments controlling Temperature decrease
		##describe<< 
		TFix=numeric(0)				##<< numeric vector (nResComp) specifying a finite value for components with fixed Temperatue, and a non-finite for others
		, TMax=numeric(0)			##<< numeric vector (nResComp) specifying a maximum temperature for result components.
		, TDecProp=0.9				##<< proportion of Temperature decrease: below one to diminish risk of decreasing Temperature too fast (below what is supported by other data streams)
		##end<<
	)
	, ctrlConvergence = list(		##<< list or arguments controlling check for convergence
		##describe<< 
		maxRelTChangeCrit=0.05 	    ##<< if Temperature of the components changes less than specified value, the algorithm can finish
        , minTBase0 = 1e-4          ##<< if Temperature-1  gets below this temperature, the algorithm can finish
		, maxLogDenDriftCrit=0.3	##<< if difference between mean logDensity of first and fourth quartile of the sample is less than this value, we do not need further batches because of drift in logDensity
		, gelmanCrit=1.4			##<< do not change Temperature, if variance between chains is too high, i.e. Gelman Diag is above this value
		, critSpecVarRatio=20		##<< if proprotion of spectral Density to Variation is higher than this value, signal problems and resort to subspaces
		, dumpfileBasename=NULL		##<< scalar string: filename to dump stack before stopping. May set to "recover"
        , maxThin=64                ##<< scalar positive integer: maximum thinning interval. If not 0 then thinning is increased on too high spectral density 
		##end<<
	)
){
	# save calling arguments to allow an continuing an interrupted run
	argsF <- as.list(sys.call())[-1]	# do not store the file name and the first two arguments
	argsF <- argsF[ !(names(argsF) %in% c("","nGen","mc")) ]  # remove positional arguments and arguments mc and nGen
	argsFEval <- lapply( argsF, eval.parent )		# remember values instead of language objects, which might not be there on a repeated call
    mc$args <- argsFEval
    finishEarlyMsg <- NULL
    #
	#-- fill in default argument values
	if( is.null(controlTwDEMC$thin) ){ 
        controlTwDEMC$thin <- 4
    }
	thin <- controlTwDEMC$thin
	frm <- formals()
	ctrlConvergence <- if( hasArg(ctrlConvergence) ) twMergeLists( eval(frm[["ctrlConvergence"]]), ctrlConvergence ) else ctrlConvergence
	ctrlBatch <- if( hasArg(ctrlBatch) ) twMergeLists( eval(frm[["ctrlBatch"]]), ctrlBatch ) else ctrlBatch	
	ctrlT <- if( hasArg(ctrlT) ) twMergeLists( eval(frm[["ctrlT"]]), ctrlT ) else ctrlT
	if( debugSequential ){
		if( 0==length(ctrlConvergence$dumpfileBasename) ) ctrlConvergence$dumpfileBasename <- "recover"
		#if( 0==length(ctrlSubspaces$argsFSplit$debugSequential) ) ctrlSubspaces$argsFSplit$debugSequential <- TRUE 
	}
	#
	nResComp <- ncol(mc$pops[[1]]$resLogDen)
	#print("saCont: before completing temperature settings."); recover()
	ctrlT$TFix <- completeResCompVec( ctrlT$TFix, colnames(mc$pops[[1]]$resLogDen) )
	iFixTemp <- which( is.finite(ctrlT$TFix) )
    iNonFixTemp <- which( !is.finite(ctrlT$TFix) )
	ctrlT$TMax <- completeResCompVec( ctrlT$TMax, colnames(mc$pops[[1]]$resLogDen) )
	iMaxTemp <- which( is.finite(ctrlT$TMax) )	
	#
	res <- mc
	.tmp.f <- function(){
		mc0 <- concatPops(resEnd)
		#mc0 <- concatPops(res)
		matplot( mc0$temp, type="l" )
		matplot( mc0$pAccept[,1,], type="l" )
		matplot( mc0$pAccept[,2,], type="l" )
		plot( as.mcmc.list(mc0), smooth=FALSE )
		tmp <- calcTemperatedLogDen(res$pops[[1]]$resLogDen[,,1], getCurrentTemp(res) )
		matplot( tmp, type="l" )
		plot(rowSums(tmp))
		#
		matplot( mc0$resLogDen[,3,], type="l" )
		matplot( mc0$parms[,"a",], type="l" )
		matplot( mc0$parms[1:20,"b",], type="l" )
		bo <- 1:10; iPop=1
		plot( mc0$parms[bo,"a",iPop], mc0$parms[bo,"b",iPop], col=rainbow(100)[twMisc::twRescale(mc0$resLogDen[bo,"parmsSparse",iPop],c(10,100))] )
	}
    #
	# sum nObs within density - assume nObs positions correspond to resLogDenComp positions
	iDens <- seq_along(mc$dInfos)
	#iDen=1
	iCompsNonFixDen <- lapply( iDens, function(iDen){ 
			irc <-  mc$dInfos[[iDen]]$resCompPos
			irc <- irc[ !(irc %in% iFixTemp) ]
		})
	nObsDen <- sapply( iDens, function(iDen){ sum( nObs[iCompsNonFixDen[[iDen]] ]) })
	TCurr <- getCurrentTemp(mc)
    #
	iBatch=1
    nGenDone <- 0       # tracking the generations already done
	#nBatch <- ceiling( nGen/(ctrlBatch$nSamplesBatch*controlTwDEMC$thin) )
	while( nGenDone < nGen){
		if((0 < length(restartFilename)) && is.character(restartFilename) && restartFilename!=""){
			resRestart.twDEMCSA = res #avoid variable confusion on load by giving a longer name
			resRestart.twDEMCSA$nGenDone <- nGenDone	# also store updated calculation of burnin time
			resRestart.twDEMCSA$args <- argsFEval	# also store updated calculation of burnin time
			save(resRestart.twDEMCSA, file=restartFilename)
			cat(paste("Saved resRestart.twDEMCSA to ",restartFilename,"\n",sep=""))
		}
        # do the diagnostics on the last halv of the chain, but continue 
		resEnd <- thin(res, start=ceiling(getNGen(res)) * ctrlBatch$pThinPast )	# neglect the first half part
        resSparse <- squeeze(res, length.out=ceiling(getNSamples(res)*ctrlBatch$pThinPast) ) # better thin the entire chain so that distant path is sparsely kept, when continuing
        #----- check for high autocorrelation wihin spaces
		mcEnd <- concatPops(stackChainsPop(resEnd))
		#plot( as.mcmc.list(mcEnd), smooth=FALSE )
		specVarRatio <- apply( mcEnd$parms, 3, function(aSamplePurePop){
				spec <- spectrum0.ar(aSamplePurePop)$spec
				varSample <- apply(aSamplePurePop, 2, var)
				spec/varSample
			})	
		ssc <- stackChainsPop(resEnd)	# combine all chains of one population
		mcl <- as.mcmc.list(ssc)
		#plot( as.mcmc.list(mcl), smooth=FALSE )
		#plot( res$pops[[1]]$parms[,"b",1] ~ res$pops[[1]]$parms[,"a",1], col=rainbow(255)[round(twRescale(-res$pops[[1]]$resLogDen[,"logDen1",1],c(10,200)))] )
		gelmanDiagRes <- try( {tmp<-gelman.diag(mcl); if(length(tmp$mpsrf)) tmp$mpsrf else tmp$psrf[1]} )	# cholesky decomposition may throw errors
		#gelmanDiagRes <- try( gelman.diag(mcl)$mpsrf )	# cholesky decomposition may throw errors
        # compare intra-chain to inter-chain gelman diag
        #iPop <- 4
        gelmanDiagPops <- sapply( 1:getNPops(resEnd), function(iPop){
             mcl <- as.mcmc.list( subPops(resEnd,iPop) )
             tmp <- try(gelman.diag(mcl), silent=TRUE)
             if( inherits(tmp,"try-error")) 999 else if(length(tmp$mpsrf)) tmp$mpsrf else tmp$psrf[1]
         })
        maxGelmanDiagPops <- max(gelmanDiagPops)
        T0 <- calcBaseTemp(TCurr, nObs, TFix=ctrlT$TFix, iNonFixTemp=iNonFixTemp)-1
        # calculate average resLogDen of last part
        resLogDen <- stackChains(concatPops(resEnd)$resLogDen)
        .resDrift <- isLogDenDrift(resLogDen, resEnd$dInfos, maxDrift=ctrlConvergence$maxLogDenDriftCrit)
        lDenLastPart <- structure( attr(.resDrift,"resTTest")[2,], names=names(resEnd$dInfos) )
        ##details<< \describe{\item{diagnostics}{
        ## in resulting twDEMC an entry diagnostics is returned. Its a list with diagnostics after each batch.
        ## Its components have been calculated on the last halv of the chains.
        ##describe<<
        newDiags <- c(
                T0=T0   ##<< Temperature (variance inflation factor) scaled for one observation minus one
                ,maxGelmanDiagPops=maxGelmanDiagPops    ##<< maximum of within population gelman diagnostics
                ,gelmanDiag=gelmanDiagRes               ##<< gelman diag calculated between populations
                ,maxSpecVarRatio=max(specVarRatio)      ##<< ratio of spectral density to sample Variance (max over all parameters). This is a measure that increases with autocorrelation.
                ,lDenLastPart
                #,attr(.resDrift,"logDenComp")
        )
        # if number of densities changed, reset columns in diagnostics (else rbind fails)
        if( (length(resSparse$diagnostics) > 0) && (length(newDiags) != ncol(resSparse$diagnostics)) ){
            resSparse$diagnostics <- cbind( resSparse$diagnostics[,1:4]
                , matrix(NA, nrow=nrow(resSparse$diagnostics), ncol=length(lDenLastPart)
                    , dimnames=list(NULL, names(lDenLastPart))) )
        }
        resSparse$diagnostics <- rbind( resSparse$diagnostics, newDiags)
        rownames( resSparse$diagnostics ) <- NULL
        ## Last components represent the average logDensity of the last 8th of the chains.
        ##details<< }}
        #
        ##details<< \describe{\item{Adaptive thinning interval (\code{maxThin})}{ 
        ##  If spectral variance is larger than its critical ratio \code{ctrlConvergence$critSpecVarRatio}
        ##  then thinning interval needs to be increased.
        ##  If doubling the thinning interval is not larger than \code{ctrlConvergence$maxThin} then a new batch 
        ##  with higher thinning is attempted. Else the SADEMC quits with reporting an error with return values component \code{failureMsg}. 
        ##}}
        if( any(specVarRatio > ctrlConvergence$critSpecVarRatio) ){
            newThin <- resSparse$thin * 2
            if( length(ctrlConvergence$maxThin) && (ctrlConvergence$maxThin > 0) && (newThin <= ctrlConvergence$maxThin) ){
                # increase thinning interval instead of reporting error
                print(paste("twDEMCSACont: too much autocorrelation (specVarRatio=",signif(max(specVarRatio),3),"). Increasing thinning interval to ",newThin,sep=""))
                resSparse <- thin(resSparse, newThin )
                controlTwDEMC$thin <- newThin
            }else{
                res <- resSparse
                res$failureMsg <- paste("twDEMCSACont: too much autocorrelation (specVarRatio=",signif(max(specVarRatio),3),"). Try using divideTwDEMC or higher thinning interval",sep=")")
                print(res$failureMsg)
                break
            }
        }
        # check if can decrease temperature, must calc TEnd and TEnd0
        if(  !inherits(gelmanDiagRes,"try-error") && 
                maxGelmanDiagPops <= ctrlConvergence$gelmanCrit &&          # each pop converged
                (gelmanDiagRes <= 1.2*maxGelmanDiagPops)                    # all pops explore the same optimum
        ){
            if( length(ctrlT$TEndFixed)==0 ){
                #
                # print("twDEMCSACont: before Temperature calculation."); recover()    
                #debug(getBestModelIndices)
                resLogDenT <- calcTemperatedLogDen(resLogDen, TCurr)
                #iBest <- getBestModelIndices(resLogDenT, resEnd$dInfos, 0.01)
                iBest <- getBestModelIndices(resLogDenT, resEnd$dInfos, ctrlT$qBest)
                SkBest <- -2*resLogDen[iBest,,drop=FALSE]
                # sapply( obs, function(obsk){ sum(obsk$sd^2) })
                #trace(calcBaseTempSk, recover) #untrace(calcBaseTempSk)
                TEnd0Target <- calcBaseTempSk(SkBest,nObs=nObs, TFix=ctrlT$TFix, iNonFixTemp=iNonFixTemp, isVerbose=ctrlT$isVerbose) -1
                TEnd <- calcStreamTemp(TEnd0Target+1, nObs, TFix=ctrlT$TFix)
    			TEnd[iMaxTemp] <- pmin(TEnd[iMaxTemp], ctrlT$TMax[iMaxTemp])	# do not increase T above TMax
    			#relTChange <- abs(TEnd - TCurr)/TEnd
    			# slower TDecrease to avoid Temperatues that are not supported by other datastreams
                TEnd0 <- T0 - ctrlT$TDecProp*(T0-TEnd0Target) 
    			TEnd <- TCurr - ctrlT$TDecProp*(TCurr-TEnd)
            } else { # ctrlT$TEndFixed is specified
                # ultimately decrease temperature to 1
                TEnd0Target <- ctrlT$TEndFixed-1
                TEnd0 <- T0 - ctrlT$TDecProp*(T0-(ctrlT$TEndFixed-1))
                TEnd <- TCurr <- calcStreamTemp(TEnd0+1, nObs, TFix=ctrlT$TFix) 
            }
            #relTChange <- (TEnd0Target - T0)/(max(ctrlConvergence$minTBase0,TEnd0))
            relTChange <- (TEnd0 - T0)/(max(ctrlConvergence$minTBase0,TEnd0))
            cat("relTChange=",signif(relTChange,2),"\n")            
            #if( (max(relTChange) <= maxRelTChange) ) 
            # recover()
            #trace(isLogDenDrift, recover )
            # may finish early, but do at least one batch
            if( relTChange > 0 ){
                ctrlT$TDecProp <- ctrlT$TDecProp * 2/3
                cat("decreasing ctrlT$TDecProp to ",signif(ctrlT$TDecProp*100,2),"%\n")
            } 
            if( (nGenDone != 0) &&
                    ((max(abs(relTChange)) <= ctrlConvergence$maxRelTChangeCrit) || T0 <= ctrlConvergence$minTBase0 ) && 
                    !.resDrift 
                    ){
                #res <- resSparse
                # keep original res
                finishEarlyMsg <- paste("twDEMCSA: Maximum Temperture change only ",signif(max(relTChange)*100,2)
                        ,"% and no drift in logDensity (lden=",paste(signif(lDenLastPart,3),collapse=","),"). Finishing early.",sep="") 
                print( finishEarlyMsg)
                print(paste("gelmanDiagPops=",signif(gelmanDiagRes,2)
                                ," max(gelmanDiagPop)=",signif(maxGelmanDiagPops,2)
                                ," max(specVarRatio)=",signif(max(specVarRatio),2)
                                ," T0=",signif(T0,2)
                                ," logDen=",paste(signif(lDenLastPart,3),collapse=",")
                                #," T=",paste(signif(TCurr,2),collapse=",")
                                , sep=""))
                break
            }
        }else{ # cannot decrease temperature
            TEnd0 <- T0
            TEnd <-TCurr
		}
		#TMin <- pmin(TMin, TEnd)
        print(paste("gelmanDiagPops=",signif(gelmanDiagRes,2)
                ," max(gelmanDiagPop)=",signif(maxGelmanDiagPops,2)
                ," max(specVarRatio)=",signif(max(specVarRatio),2)
                ," T0=",signif(T0,2)
                ," logDen=",paste(signif(lDenLastPart,3),collapse=",")
                #," T=",paste(signif(TCurr,2),collapse=",")
            , sep=""))
        .nGenBatch <- min( ctrlBatch$nSamplesBatch*controlTwDEMC$thin,  nGen - nGenDone )
        print(paste("doing ",.nGenBatch," generations to TEnd0=",signif(TEnd0,2)
                            ,", TEnd=",paste(signif(TEnd,2),collapse=","),sep="") )
        res1 <- res <- twDEMCBlock( resSparse
			, nGen=.nGenBatch
			, debugSequential=debugSequential
			, m0=m0
			, controlTwDEMC=controlTwDEMC
			, TEnd=TEnd
			# replace lines below later on by ...
			#, blocks = blocks
			,...
		)
		TCurr <- getCurrentTemp(res)
        nGenDone <- nGenDone + .nGenBatch
		print(paste("finished ",nGenDone," out of ",nGen," gens.   ", date(), sep="") )
	} # end while( nGenDone < nGen)
	res$TGlobal <- max(getCurrentTemp(res)[unlist(iCompsNonFixDen)])
	res$args <- argsFEval
    res$finishEarly <- finishEarlyMsg
    #----------------------  update diagnostics for last batch (replicated from start of the loop)
    resEnd <- thin(res, start=ceiling(getNGen(res)) * ctrlBatch$pThinPast )	# neglect the first half part
    mcEnd <- concatPops(stackChainsPop(resEnd))
    #----- check for high autocorrelation wihin spaces
    #plot( as.mcmc.list(mcEnd), smooth=FALSE )
    specVarRatio <- apply( mcEnd$parms, 3, function(aSamplePurePop){
                spec <- spectrum0.ar(aSamplePurePop)$spec
                varSample <- apply(aSamplePurePop, 2, var)
                spec/varSample
            })	
    ssc <- stackChainsPop(resEnd)	# combine all chains of one population
    mcl <- as.mcmc.list(ssc)
    #plot( as.mcmc.list(mcl), smooth=FALSE )
    #plot( res$pops[[1]]$parms[,"b",1] ~ res$pops[[1]]$parms[,"a",1], col=rainbow(255)[round(twRescale(-res$pops[[1]]$resLogDen[,"logDen1",1],c(10,200)))] )
    gelmanDiagRes <- try( {tmp<-gelman.diag(mcl); if(length(tmp$mpsrf)) tmp$mpsrf else tmp$psrf[1]} )	# cholesky decomposition may throw errors
    #gelmanDiagRes <- try( gelman.diag(mcl)$mpsrf )	# cholesky decomposition may throw errors
    # compare intra-chain to inter-chain gelman diag
    #iPop <- 4
    gelmanDiagPops <- sapply( 1:getNPops(resEnd), function(iPop){
                mcl <- as.mcmc.list( subPops(resEnd,iPop) )
                tmp <- try(gelman.diag(mcl), silent=TRUE)
                if( inherits(tmp,"try-error")) 999 else if(length(tmp$mpsrf)) tmp$mpsrf else tmp$psrf[1]
            })
    maxGelmanDiagPops <- max(gelmanDiagPops)
    T0 <- calcBaseTemp(TCurr, nObs, TFix=ctrlT$TFix, iNonFixTemp=iNonFixTemp)-1
    # calculate average resLogDen of last part
    resLogDen <- stackChains(concatPops(resEnd)$resLogDen)
    .resDrift <- isLogDenDrift(resLogDen, resEnd$dInfos, maxDrift=ctrlConvergence$maxLogDenDriftCrit)
    lDenLastPart <- structure( attr(.resDrift,"resTTest")[2,], names=names(resEnd$dInfos) )
    newDiags <- c(
            T0=T0   ##<< Temperature (variance inflation factor) scaled for one observation minus one
            ,maxGelmanDiagPops=maxGelmanDiagPops    ##<< maximum of within population gelman diagnostics
            ,gelmanDiag=gelmanDiagRes               ##<< gelman diag calculated between populations
            ,maxSpecVarRatio=max(specVarRatio)      ##<< ratio of spectral density to sample Variance (max over all parameters). This is a measure that increases with autocorrelation.
            ,lDenLastPart
    #,attr(.resDrift,"logDenComp")
    )
    res$diagnostics <- rbind( res$diagnostics, newDiags)
    rownames( res$diagnostics ) <- NULL
	res
}

.tmp.f <- function(){       # ggplotResults
	mc0 <- concatPops(res)
	.nSample <- 128
	dfDen <- rbind(
		cbind( data.frame( scenario="S1", {tmp <- stackChains(mc0)[,-(1:getNBlocks(mc0))]; tmp[ seq(1,nrow(tmp),length.out=.nSample),] }))
	)
	dfDenM <- melt(dfDen)
	#str(dfDenM)
	require(ggplot2)
	StatDensityNonZero <- StatDensity$proto( calculate<- function(.,data, scales, dMin=0.001, ...){
			res <- StatDensity$calculate(data,scales,...)
			# append a zero density at the edges
			if( res$density[1] > 0) res <- rbind( data.frame(x=res$x[1]-diff(res$x[1:2]),density=0,scaled=0,count=0),res)
			if( res$density[nrow(res)] > 0) res <- rbind( res, data.frame(x=res$x[nrow(res)]+diff(res$x[nrow(res)-c(1,0)]),density=0,scaled=0,count=0))
			bo <- (res$density < dMin) & (c(res$density[-1],0) < dMin) & (c(0,res$density[-length(res$density)]) < dMin)
			res$density[ bo ] <- NA
			res$count[ bo ] <- NA
			res$scaled[ bo ] <- NA
			res
		}, required_aes=c("x"))
#with(StatDensityNonZero, mtrace(calculate)) 
#with(StatDensityNonZero, mtrace(calculate,F)) 
	
	stat_densityNonZero <- function (
		### constructs a new StatPrior statistics based on aesthetics x and parName
		mapping = NULL, data = NULL, geom = "line", position = "stack", 
		adjust = 1,	...
	){ 
		StatDensityNonZero$new(mapping = mapping, data = data, geom = geom, 
			position = position, adjust = adjust, ...)
	}
	
#pgm <- geom_ribbon( alpha=0.8, aes(ymax = ..density.., ymin = -..density..), stat = "density")
#pgm <- geom_ribbon(  aes(ymax = ..density.., ymin = -..density..), stat = "density")
#pgm <- stat_densityNonZero(  aes(ymax = -..density..),  size=1, geom="line")
	pgm <- geom_line(aes(y = +..scaled..), stat="densityNonZero", size=1)
	pgm2 <- geom_line(aes(y = -..scaled..), stat="densityNonZero", size=1) 
	optsm <- opts(axis.title.x = theme_blank(), axis.text.y = theme_blank(), axis.ticks.y = theme_blank() ) 
	pa <- ggplot(dfDen, aes(x = a, colour=scenario, linetype=scenario)) + pgm + pgm2 + optsm + scale_y_continuous('a')
	pb <- ggplot(dfDen, aes(x = b, colour=scenario, linetype=scenario)) + pgm + pgm2 + optsm + scale_y_continuous('b')
#pb
	windows(width=7, height=3)
	grid.newpage()
	pushViewport( viewport(layout=grid.layout(2,1)))
	print(pa , vp = viewport(layout.pos.row=1,layout.pos.col=1))	
	print(pb + opts(legend.position = "none") , vp = viewport(layout.pos.row=2,layout.pos.col=1))	
	
	#----------- ggplot predictive posterior
	#scenarios <- c("R","RS","RSw","DG","DM")
	scenarios <- c("R")
	nScen <- length(scenarios)
	# infer quantiles of predictions
	.nSample=128
	y1M <- array( NA_real_, dim=c(.nSample, length(twTwoDenEx1$obs$y1),nScen), dimnames=list(sample=NULL,iObs=NULL,scenario=scenarios)  )
	y2M <- array( NA_real_, dim=c(.nSample, length(twTwoDenEx1$obs$y2),nScen), dimnames=list(sample=NULL,iObs=NULL,scenario=scenarios) )
	resScen <- #list( res1, res2, res2b, res3a, res3 ); names(resScen) <- scenarios
	resScen <- list( mc0 ); names(resScen) <- scenarios
	#scen <- "R"
	for( scen in scenarios ){
		ss <- stackChains(resScen[[scen]])[,-(1:getNBlocks(resScen[[scen]]))]
		ssThin <- ss[round(seq(1,nrow(ss),length.out=.nSample)),]
		#i <- .nSample
		for( i in 1:.nSample){
			pred <-  twTwoDenEx1$fModel(ssThin[i,], xSparse=twTwoDenEx1$xSparse, xRich=twTwoDenEx1$xRich, thresholdCovar=thresholdCovar) 
			y1M[i,,scen] <- pred$y1
			y2M[i,,scen] <- pred$y2
		}
	}
	predQuantilesY1 <- lapply( scenarios, function(scen){
			t(apply( y1M[,,scen],2, quantile, probs=c(0.025,0.5,0.975) ))
		}); names(predQuantilesY1) <- scenarios
	predQuantilesY2 <- lapply( scenarios, function(scen){
			t(apply( y2M[,,scen],2, quantile, probs=c(0.025,0.5,0.975) ))
		}); names(predQuantilesY2) <- scenarios
	str(predQuantilesY1)
	rm( y1M, y2M)	# free space, we only need the summary
	
	#str(pred1)
	dfPred <- rbind(
		cbind( data.frame( scenario="R", variable="y1", value = pred1$y1 ), predQuantilesY1[["R"]] )
		,cbind( data.frame( scenario="R", variable="y2", value = pred1$y2 ), predQuantilesY2[["R"]] )
	)
	dfPred$observations <- NA_real_
	dfPred$observations[dfPred$variable=="y1"] <- twTwoDenEx1$obs$y1
	dfPred$observations[dfPred$variable=="y2"] <- twTwoDenEx1$obs$y2
	colnames(dfPred)[ match(c("2.5%","50%","97.5%"),colnames(dfPred)) ] <- c("lower","median","upper")
	str(dfPred)
	
	p1 <- ggplot(dfPred, aes(x=value, y=observations, colour=scenario) ) +
		geom_errorbarh(aes(xmax = upper, xmin = lower)) +
		geom_point() +
		facet_wrap( ~ variable, scales="free") +
		opts(axis.title.x=theme_blank() ) +
		geom_abline(colour="black") +
		c()
	p1
	
	
}

.tmp.f <- function(){ # does not work        # AllDensities
	# same as example but with each parameter updated against both densities
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for updating each in turn against both densities
	dInfos=list(
		dSparse=list(fLogDen=denSparsePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		bSparse=list(dInfoPos="dSparse", compPos=c("a","b") )
		,bRich=list(dInfoPos="dRich", compPos=c("a","b") )
	)
	
	do.call( dInfos$dSparse$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparse$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	
	#str(twTwoDenEx1)
	nObs <- c( obsSparse=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	#trace(twDEMCSA, recover )
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs, debugSequential=TRUE, ctrlBatch=list(nSamplesBatch=64) )
	# continue run for 512 generations
	res <- twDEMCSA( res, 512 )
	
	(TCurr <- getCurrentTemp(resPops))
	mc0 <- concatPops(res)
	logDenT <- calcTemperatedLogDen(stackChains(mc0$resLogDen), TCurr)
	iBest <- getBestModelIndex( logDenT, res$dInfos )
	maxLogDenT <- logDenT[iBest, ]
	ss <- stackChains(mc0$parms)
	(thetaBest <- ss[iBest, ])
	(.qq <- apply(ss,2,quantile, probs=c(0.025,0.5,0.975) ))
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(rowSums(logDenT),c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"obsSparse"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"logDen1"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rgb(
			twRescale(logDenT[,"obsSparse"]),0, twRescale(logDenT[,"logDen1"]) ))
	apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )
	
	# density of parameters
	plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")
	
	# predictive posterior (best model only)
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparse=xSparse, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
	
}

.tmp.f <- function(){ #         # AllSparse
	# same as example but with b also updated against Sparse, a only against Sparse 
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for updating each in turn against both densities
	dInfos=list(
		dSparse=list(fLogDen=denSparsePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		bSparse=list(dInfoPos="dSparse", compPos=c("a","b") )
		,bRich=list(dInfoPos="dRich", compPos=c("a") )
	)
	
	do.call( dInfos$dSparse$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparse$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	#str(twTwoDenEx1)
	nObs <- c( obsSparse=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	#untrace(twDEMCSA )
	#trace(twDEMCSA, recover )
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs, nBatch=5, debugSequential=TRUE, TFix=c(1,NA,NA) )
	resPops <- res <- twDEMCSACont( res, 256, nObs=nObs, debugSequential=TRUE
		, ctrlT=list( TFix=c(parmsSparse=1) ) 
	)
	
	(TCurr <- getCurrentTemp(resPops))
	mc0 <- concatPops(res)
	logDenT <- calcTemperatedLogDen(stackChains(mc0$resLogDen), TCurr)
	iBest <- getBestModelIndex( logDenT, res$dInfos )
	maxLogDenT <- logDenT[iBest, ]
	ss <- stackChains(mc0$parms)
	(thetaBest <- ss[iBest, ])
	(.qq <- apply(ss,2,quantile, probs=c(0.025,0.5,0.975) ))
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(rowSums(logDenT),c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"obsSparse"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"logDen1"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rgb(
			twRescale(logDenT[,"obsSparse"]),0, twRescale(logDenT[,"logDen1"]) ))
	apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )
	
	# density of parameters
	plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")
	
	# predictive posterior (best model only)
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparse=xSparse, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
	
}


.tmp.f <- function(){ # does not work        # .tmp.2DenSwitched
	# same as example but with each parameter b updated against Sparse
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for updating each in turn against both densities
	dInfos=list(
		dSparse=list(fLogDen=denSparsePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		bSparse=list(dInfoPos="dSparse", compPos=c("b") )
		,bRich=list(dInfoPos="dRich", compPos=c("a") )
	)
	
	do.call( dInfos$dSparse$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparse$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	
	#str(twTwoDenEx1)
	nObs <- c( obsSparse=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	#untrace(twDEMCSA )
	#trace(twDEMCSA, recover )
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs, nBatch=5, debugSequential=TRUE )
	
	(TCurr <- getCurrentTemp(resPops))
	mc0 <- concatPops(res)
	logDenT <- calcTemperatedLogDen(stackChains(mc0$resLogDen), TCurr)
	iBest <- getBestModelIndex( logDenT, res$dInfos )
	maxLogDenT <- logDenT[iBest, ]
	ss <- stackChains(mc0$parms)
	(thetaBest <- ss[iBest, ])
	(.qq <- apply(ss,2,quantile, probs=c(0.025,0.5,0.975) ))
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(rowSums(logDenT),c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"obsSparse"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"logDen1"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rgb(
			twRescale(logDenT[,"obsSparse"]),0, twRescale(logDenT[,"logDen1"]) ))
	apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )
	
	# density of parameters
	plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")
	
	# predictive posterior (best model only)
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparse=xSparse, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
	
}

.tmp.f <- function(){   # oneDensity_andCluster
	# same as example but with only one combined density for both parameters
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for only one logDensity - Temperature of the strongest component goes to zero and other increase according to mismatch
	dInfos=list( list(fLogDen=denBoth, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta)) )
	
	do.call( dInfos[[1]]$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos[[1]]$argsFLogDen))
	#str(twTwoDenEx1)
	nObs <- c( y1=length(twTwoDenEx1$obs$y1), y2=length(twTwoDenEx1$obs$y2) )
	
	#trace(twDEMCSA, recover)
	argsTwDEMCSA <- list( thetaPrior=thetaPrior, covarTheta=covarTheta, dInfos=dInfos, nObs=nObs
		,TFix = c(1,NA,NA)
		,nGen=256
		#,TMax = c(NA,1.2,NA)	# do not increase Sparse observations again too much
		, nBatch=3 
		, debugSequential=TRUE
	)
	resPops <- res <- do.call( twDEMCSA, argsTwDEMCSA )
	
	.tmp.f <- function(){ #continueRun
		res$pops[[2]]$spaceInd <- 4		# test spaces different from 1:nSpace
		res2 <- twDEMCSA(res)
	}
	.tmp.f <- function(){ #byCluster
		runClusterParms <- list(
			fSetupCluster = function(){library(twDEMC)}
			,fRun = twDEMCSA
			,argsFRun = within( argsTwDEMCSA, debugSequential<-FALSE )
		)
		save(runClusterParms, file=file.path("..","..","projects","asom","parms","saOneDensity.RData"))
		# from asom directory run: 
		#   R CMD BATCH --vanilla '--args paramFile="parms/saOneDensity.RData"' runCluster.R Rout/runCluster_0.rout 
		# or from libra home directory run 
		#   ./bsubr_i.sh runCluster.R iproc=0 nprocSinge=1 'paramFile="parms/saOneDensity.RData"'
		# or
		#    bsub -q SLES -n 4 ./bsubr_i.sh runCluster.R  nprocSingle=4 'paramFile="parms/saOneDensity.RData"'
		load("parms/res_saOneDensity_1.RData")
	}
	
	(TCurr <- getCurrentTemp(res))
	mc0 <- concatPops(res)
	logDenT <- calcTemperatedLogDen(stackChains(mc0$resLogDen), TCurr)
	iBest <- getBestModelIndex( logDenT, res$dInfos )
	maxLogDenT <- logDenT[iBest, ]
	ss <- stackChains(mc0$parms)
	(thetaBest <- ss[iBest, ])
	(.qq <- apply(ss,2,quantile, probs=c(0.025,0.5,0.975) ))
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(rowSums(logDenT),c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"y1"],c(10,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"y2"],c(10,100))] )
	plot( ss[,"a"], ss[,"b"], col=rgb(
			twRescale(logDenT[,"y1"]),0, twRescale(logDenT[,"y2"]) ))
	apply( apply( logDenT, 2, quantile, probs=c(0.2,0.8) ),2, diff )
	
	# density of parameters
	plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")
	
	# predictive posterior (best model only)
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparse=xSparse, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
}




isLogDenDrift <- function(
	### check whether first quartile all the logDensities is significantly smaller than last quartile 
	logDenT		##<< numeric array (nStep x nResComp): logDensity (highest are best)
	, dInfos 	##<< list of lists with entry resCompPos (integer vector) specifying the position of result components for each density
	, alpha=0.05	##<< the significance level for a difference
	, maxDrift=0.3	##<< difference in LogDensity, below which no drift is signalled
){
	##details<<
	## Because of large sample sizes, very small differences may be significantly different.
	## Use argument minDiff to specify below which difference a significant difference is not regarded as drift.
	nr <- nrow(logDenT)
	nr4 <- nr %/% 4
    ds1 <- logDenT[1:nr4,   ,drop=FALSE]        # first quater of the dataset
    ds4 <- logDenT[(nr+1-nr4):nr,  ,drop=FALSE]     # last quater of the dataset
	#iDen <- 1
	tres <- sapply( seq_along(dInfos), function(iDen){
			dInfo <- dInfos[[iDen]]
			rs1 <- rowSums(ds1[,dInfo$resCompPos ,drop=FALSE])
			rs4 <- rowSums(ds4[,dInfo$resCompPos ,drop=FALSE])
			resTTest <- try(t.test(rs1,rs4,"less"), silent=TRUE)
            if( inherits(resTTest,"try-error") ){  # logDen essentially constant
    			#bo <- (diff(resTTest$estimate) >= maxDrift ) && (resTTest$p.value <= alpha)
                c(rs1,rs4,1)
            } else {
                c( resTTest$estimate, p=resTTest$p.value)
            }   
		})
    #tresi <- tres[,1]
    boL <- apply( tres, 2, function(tresi){  
                (diff(tresi[1:2]) >= maxDrift ) && (tresi["p"] <= alpha)
            })
	ret <- any( boL )
    attr( ret, "resTTest") <- tres
    attr( ret, "logDenComp") <- colMeans(ds4)       # average logDensity of all components of last quarter
    ##value<< TRUE if any of the logDensities are a significantly greater in the fourth quantile compared to the first quantile of the samples
    ## , attribute \code{resTTest}: numeric maxtrix (logDenStart, logDenEnd, pTTest x nDInfo )
    ## , attribute \code{logDenComp}: numeric vector (nResComp): average logDenT of all components of last quater
    ret
}
attr(isLogDenDrift,"ex") <- function(){
	data(twdemcEx1)
	logDenT <- calcTemperatedLogDen( twdemcEx1, getCurrentTemp(twdemcEx1))
    #undebug(isLogDenDrift)
	isLogDenDrift(logDenT, twdemcEx1$dInfos )
}


twRunDEMCSA <- function(
	### wrapper to twDEMCSA compliant to runCluster.R
	argsTwDEMCSA	 
	### Arguments passed to twDEMCSA -> twDEMCBlockInt
	### It is updated by \dots.
	,...				 	##<< further arguments passed to twDEMCSA -> twDEMCBlockInt
	,prevResRunCluster=NULL	##<< results of call to twRunDEMCSA, argument required to be called from runCluster.R
	,restartFilename=NULL	##<< name of the file to store restart information, argument required to be called from runCluster.R 
){
	require(twDEMC)
	#update argsDEMC to args given 
	argsDEMC <- argsTwDEMCSA
	# update the restartFilename argument
	.dots <- list(...)
	argsDEMC[ names(.dots) ] <- .dots
	#for( argName in names(.dots) ) argsDEMC[argName] <- .dots[argName]
	if( 0 != length(restartFilename)) argsDEMC$restartFilename <- restartFilename
	# do the actual call
	res <- do.call( twDEMCSA, argsDEMC )
}
.tmp.f <- function(){       # testRestart
	rFileName <- file.path("tmp","example_twDEMCSA.RData")
	load(rFileName)	# resRestart.twDEMCSA
	#untrace(twDEMCSACont)
	#trace(twDEMCSACont, recover )
	res1 <- do.call( twDEMCSACont, c( list(mc=resRestart.twDEMCSA), resRestart.twDEMCSA$args ))
}

completeResCompVec <- function(
	### given a vector with names, make a vector corresponding to resComp with components, not yet given filled by NA
	x				##<< named vector to be completed
	,resCompNames	##<< names of the result components
){
	xAll <- x
	nResComp <- length(resCompNames)
	if( 0 == length( x) ) 
		xAll <- structure( rep( NA_real_, nResComp), names=resCompNames )
	#else if( nResComp != length( x) ){ # twutz 141120: can give a different names of same length, must correct
    else{
		iFix <- sapply( names(x), match, resCompNames )
		if( any(is.na(iFix))){
            misNames <- paste(names(x)[ is.na(iFix)], collapse=",")
            stop(paste("completeResCompVec: x must be of the same length as resCompNames or all its names must correspond names of resCompNames. Unmatched names:",misNames))            
        } 
		xAll <- structure( rep( NA_real_, nResComp), names=resCompNames )
        xfin <- x[ is.finite(x) ]
        #iName <- names(xfin)[1]
		for( iName in names(xfin) ){
            xAll[ names(xAll)== iName ] <- xfin[iName]
        }
	}
	xAll
	### vector with names resCompNames, with corresponding values of x set, others NA
}







