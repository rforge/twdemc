# slowly decreasing Temperature (cost reduction factor)

twDEMCSA <- function(
	### simulated annealing DEMC 
	thetaPrior			##<< vector of parameters, point estimate
	,covarTheta			##<< the a prior covariance of parameters
	,nGen=512			##<< number of generations in the initial batch, default 512
	,nObs				##<< integer vector (nResComp) specifying the number of observations for each result component
	,... 				##<< further argument to \code{\link{twDEMCBlockInt}}
	, dInfos			##<< argument to \code{\link{twDEMCBlockInt}}
	, m0 = calcM0twDEMC(length(thetaPrior),nChainPop)	##<< minimum number of samples in step for extending runs
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMCBlockInt}} containing entry thin
	, debugSequential=FALSE		##<< set to TRUE to avoid parallel execution, good for debugging
	, restartFilename=NULL		##<< filename to write intermediate results to
	#
	,nChainPop=4		##<< number of chains within population
	,nPop=2				##<< number of populations
	,doIncludePrior=FALSE	##<< should the prior be part of initial population 
		##<< Recommendation to set to false, because if the TRUE parameter is in initial set, the Temperature is set to 1
	#
	, ctrlBatch = list(             ##<< list of arguments controlling batch executions, see \code{\link{twDEMCSACont}} and \code{\link{divideTwDEMCSACont}}
		##describe<< 
		nGen0=m0*controlTwDEMC$thin*3	##<< number of generations for the intial run
		,useSubspaceAdaptation=FALSE	##<< if TRUE then overall space is devided and each subspace is explored with locally adapted DEMC, see	\code{\link{divideTwDEMCSACont}}	
		##end<<
		)	
	, ctrlT = list(                 ##<< list of arguments controlling Temperature decrease, see \code{\link{twDEMCSACont}}  and \code{\link{divideTwDEMCSACont}}
		TFix=numeric(0)				##<< numeric vector (nResComp) specifying a finite value for components with fixed Temperatue, and a non-finite for others
		, TMax=numeric(0)			##<< numeric vector (nResComp) specifying a maximum temperature for result components.
		##describe<< 
		, qTempInit=0.4		##<< quantile of logDensities used to calculate initial beginning and end temperature, with default 0.4: 40% of the space is accepted 
		##end<<
	)
	, ctrlConvergence = list()		##<< list or arguments controlling check for convergence, see \code{\link{twDEMCSACont}}  and \code{\link{divideTwDEMCSACont}}
	, ctrlSubspaces = list()		##<< list of arguments controlling splitting and merging of subspaces, see \code{\link{divideTwDEMCSACont}}
){
	##detail<< \describe{\item{Continuing a previous run}{
	## When supplying the first of class twDEMCPops to twDEMCSA, the run is extended.
	## All other parameters except nGen, then are ignored.
	## In order to change parameters, modify list entry args in the twDEMCPops object.
	##}}
	# in order to continue a previous result, supply it as a first argument and add entry args
	if( inherits(thetaPrior,"twDEMCPops") && 0 != length(thetaPrior$args) ){
		if( 0!=length(restartFilename) )
			thetaPrior$args$restartFilename=restartFilename
		ret <- if( 0!=length(thetaPrior$args$ctrlBatch) &&
				   0!=length(thetaPrior$args$ctrlBatch$useSubspaceAdaptation) &&	
				   thetaPrior$args$ctrlBatch$useSubspaceAdaptation 
			)
			do.call(divideTwDEMCSACont,c(list(thetaPrior, nGen=nGen),thetaPrior$args))
		else	
			do.call(twDEMCSACont,c(list(thetaPrior, nGen=nGen),thetaPrior$args))
		#ret$args$ctrlBatch$useSubspaceAdaptation <- ctrlBatch$useSubspaceAdaptation
		return(ret)
	}
	#-- fill in default argument values
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
	thin <- controlTwDEMC$thin
	frm <- formals()
	ctrlConvergence <- if( hasArg(ctrlConvergence) ) twMergeLists( eval(frm[["ctrlConvergence"]]), ctrlConvergence ) else ctrlConvergence
	ctrlBatch <- if( hasArg(ctrlBatch) ) twMergeLists( eval(frm[["ctrlBatch"]]), ctrlBatch ) else ctrlBatch	
	ctrlT <- if( hasArg(ctrlT) ) twMergeLists( eval(frm[["ctrlT"]]), ctrlT ) else ctrlT
	#
	#mtrace(initZtwDEMCNormal)
	Zinit0 <- initZtwDEMCNormal( thetaPrior, covarTheta, nChainPop=nChainPop, nPop=nPop, doIncludePrior=doIncludePrior, m0=m0)
	ss <- stackChains(Zinit0)
	logDenL <- lapply( dInfos, function(dInfo){
		resLogDen <- twCalcLogDenPar( dInfo$fLogDen, ss, argsFLogDen=dInfo$argsFLogDen)$logDenComp
	})
	#dInfo <- dInfos[[1]]
	#tmp <- do.call( dInfo$fLogDen, c(list(ss[1,]), dInfo$argsFLogDen) )
	# replace missing cases
	logDenDS0 <- abind( logDenL )
	# logDenDS0[1,1] <- NA		# testing replacement of non-finite cases
	# logDenDS0[1:round(nrow(logDenDS0)/1.8),1] <- NA		# testing replacement of non-finite cases
	boFinite <- apply(is.finite(logDenDS0), 1, all )
	m0FiniteFac <- sum(boFinite) / nrow(logDenDS0)
	if( m0FiniteFac == 1){
			Zinit <- Zinit0
			logDenDS <- logDenDS0
		}else if( m0FiniteFac < 0.05 ){
			stop("twDEMCSA: less than 5% finite solutions in initial exploration. Check setup.")
		}else if( m0FiniteFac > 0.9 ){
			sumLogDen0 <- rowSums(logDenDS0)
			Zinit <- replaceZinitNonFiniteLogDens( Zinit0, sumLogDen0 )
			logDenDS <- logDenDS0
			#XXTODO: adjust logDen by actual replacement state
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
			# logDenDS1[ 1:round(nrow(logDenDS1)/1.7),1] <- NA		# testing replacement of non-finite cases
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
	# now the legnth of resComp is known (colnames(logDenDS)), so check the TFix and TMax parameters
	if( 0==length(colnames(logDenDS)) ) 
		colnames(logDenDS) <- paste("den",1:ncol(logDenDS), sep="")
	nResComp <- ncol(logDenDS)
	ctrlT$TFix <- completeResCompVec( ctrlT$TFix, colnames(logDenDS) )
	iFixTemp <- which( is.finite(ctrlT$TFix) )
	ctrlT$TMax <- completeResCompVec( ctrlT$TMax, colnames(logDenDS) )
	iMaxTemp <- which( is.finite(ctrlT$TMax) )	
	#
	#if( 0 == length( TMaxInc) ) TMaxInc <- structure( rep( NA_real_, nResComp), names=colnames(logDenDS) )
	#if( nResComp != length( TMaxInc) ) stop("twDEMCSA: TMaxInc must be of the same length as number of result Components.")
	#iMaxIncTemp <- which( is.finite(TMaxInc) )
	#
	#------ intial temperatures: 
	#sapply( seq_along(expLogDenBest), function(i){ max(logDenDS[,i]) })
	##details<< \describe{\item{initial temperature}{
	## Initial parameters are ranked according to their maximum log-Density across components.
	## The parameters at rank positions defined by qTempInit are selected.
	## The minimum log-Density or each datastream is used obtain initial temperature as: T =  -2/nObs logDen
	##}}
	#print("twDEMCSA: before calculating initial temperature"); recover()
	rankLogDenDS <- apply(-logDenDS, 2, rank)	# highest logDen first
	iNonFixTemp <- (1:nResComp)[ ifelse(0!=length(iFixTemp),-iFixTemp, TRUE) ]
	maxRankLogDen <- apply(rankLogDenDS[ ,iNonFixTemp ,drop=FALSE],1,max)	
	# quantile produces intermediate values, which may round wrong
	# search the component nearest to the quantile (may not match exact because of duplicates)
	iQuant <- which.min( abs(maxRankLogDen - round(ctrlT$qTempInit*length(maxRankLogDen))) )
	qLogDenDS <- apply( logDenDS[iQuant, ,drop=FALSE], 2, min )
	temp0 <- tempQ <-  pmax(1,-2/nObs* qLogDenDS)
	#names(temp0) <- colnames(logDenDS)
	ctrlT$TMax[iNonFixTemp] <- pmin(ifelse(is.finite(ctrlT$TMax),ctrlT$TMax, Inf), ifelse(is.finite(tempQ),tempQ,Inf) )[iNonFixTemp]		# decrease TMax  
	temp0[iFixTemp] <- ctrlT$TFix[iFixTemp]
	if( !all(is.finite(temp0)) ) stop("twDEMCSA: encountered non-finite Temperatures.")
	#if( any(temp0 > 8)) stop("twDEMCSA: encountered too high temperature.")	
	print(paste("initial T=",paste(signif(temp0,2),collapse=","),"    ", date(), sep="") )
	#if( sum(temp0 > 1) == 0){ dump.frames("parms/debugDump",TRUE); stop("twDEMCSA: no temperature > 0") }	
	#
	ret <- res0 <-  res <- twDEMCBlock( Zinit
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
	.tmp.f <- function(){
		windows(record=TRUE)
		mc0 <- concatPops(res)
		plot( as.mcmc.list(mc0), smooth=FALSE )
		matplot( mc0$temp, type="l" )
		matplot( mc0$pAccept[,1,], type="l" )
		matplot( mc0$pAccept[,2,], type="l" )
		#trace(calcTemperatedLogDen.default, recover)
		tmp <- calcTemperatedLogDen(adrop(res$pops[[1]]$resLogDen[,,1 ,drop=FALSE],3), getCurrentTemp(res) )
		matplot( tmp, type="l" )
		plot(rowSums(tmp))
		#
		matplot( mc0$resLogDen[,3,], type="l" )
		matplot( mc0$parms[,"a",], type="l" )
		matplot( mc0$parms[1:20,"b",], type="l" )
		bo <- 1:10; iPop=1
		plot( mc0$parms[bo,"a",iPop], mc0$parms[bo,"b",iPop], col=rainbow(100)[twRescale(mc0$resLogDen[bo,"parmsSparce",iPop],c(10,100))] )
	}
	print(paste("finished initial ",ctrlBatch$nGen0," out of ",nGen," gens.    ", date(), sep="") )
	#
	if( nGen > ctrlBatch$nGen0 ){
		#nObsLocal <- nObs 
		ret <- if( ctrlBatch$useSubspaceAdaptation ){
			divideTwDEMCSACont( mc=res0, nGen=nGen-ctrlBatch$nGen0, nObs=nObs
				, m0=m0, controlTwDEMC=controlTwDEMC, debugSequential=debugSequential, restartFilename=restartFilename
				, ctrlBatch=ctrlBatch, ctrlT=ctrlT, ctrlConvergence=ctrlConvergence, ctrlSubspaces=ctrlSubspaces
				, ...
			) 
		}else {
			twDEMCSACont( mc=res0, nGen=nGen-ctrlBatch$nGen0, nObs=nObs
				, m0=m0, controlTwDEMC=controlTwDEMC, debugSequential=debugSequential, restartFilename=restartFilename
				, ctrlBatch=ctrlBatch, ctrlT=ctrlT, ctrlConvergence=ctrlConvergence
				, ...
				)
		}
	}	
	#ret$args$ctrlBatch$useSubspaceAdaptation <- ctrlBatch$useSubspaceAdaptation	# remember this argument
	ret
}
attr(twDEMCSA,"ex") <- function(){
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	dInfos=list(
		dSparce=list(fLogDen=denSparcePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		a=list(dInfoPos="dSparce", compPos="a")
		,b=list(dInfoPos="dRich", compPos="b")
	)
	
	do.call( dInfos$dSparce$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparce$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	
	#str(twTwoDenEx1)
	nObs <- c( parmsSparce=1, obsSparce=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	#trace(twDEMCSACont, recover )
	#trace(twDEMCSA, recover )
	res <- res0 <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs
		, nGen=3*256
		, ctrlT=list( TFix=c(parmsSparce=1) )
		, ctrlBatch=list( nGenBatch=256 )
		, debugSequential=TRUE
		#, restartFilename=file.path("tmp","example_twDEMCSA.RData")
	)
	res <- twDEMCSA( res0, nGen=4*256 )	# extend the former run

	(TCurr <- getCurrentTemp(res))
	mc0 <- concatPops(res)
	plot( as.mcmc.list(mc0) , smooth=FALSE )
	matplot( mc0$temp, type="l" )
	logDenT <- calcTemperatedLogDen(stackChains(mc0$resLogDen), TCurr)
	iBest <- getBestModelIndex( logDenT, res$dInfos )
	maxLogDenT <- logDenT[iBest, ]
	ss <- stackChains(mc0$parms)
	(thetaBest <- ss[iBest, ])
	(.qq <- apply(ss,2,quantile, probs=c(0.025,0.5,0.975) ))
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(rowSums(logDenT),c(10,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"obsSparce"],c(10,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"logDen1"],c(10,100))] )
	plot( ss[,"a"], ss[,"b"], col=rgb(
			twRescale(logDenT[,"obsSparce"]),0, twRescale(logDenT[,"logDen1"]) ))
	apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )
	
	# density of parameters
	plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")
	
	# predictive posterior (best model only)
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
	
}

twDEMCSACont <- function(
	### continuing simulated annealing DEMC based on previous reslt
	mc					##<< result of twDEMCBlock
	,nGen=512			##<< number of generations in the initial batch, default 512
	,nObs				##<< integer vector (nResComp) specifying the number of observations for each result component
	,... 				##<< further argument to \code{\link{twDEMCBlockInt}}
	, m0 = calcM0twDEMC(getNParms(mc),getNChainsPop(mc))	##<< minimum number of samples in step for extending runs
	, controlTwDEMC = list()		##<< list argument to \code{\link{twDEMCBlockInt}} containing entry thin
	, debugSequential=FALSE		##<< set to TRUE to avoid parallel execution, good for debugging
	, restartFilename=NULL		##<< filename to write intermediate results to
	#
	, ctrlBatch = list(				##<< list of arguments controlling batch executions
		##describe<< 
		nGenBatch=m0*controlTwDEMC$thin*5		##<< number of generations for one call to twDEMCStep
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
		maxRelTChangeCrit=0.025 	##<< if Temperature of the components changes less than specified value, the algorithm can finish
		, maxLogDenDriftCrit=0.3	##<< if difference between mean logDensity of first and fourth quartile of the sample is less than this value, we do not need further batches because of drift in logDensity
		, gelmanCrit=1.4			##<< do not change Temperature, if variance between chains is too high, i.e. Gelman Diag is above this value
		, critSpecVarRatio=20		##<< if proprotion of spectral Density to Variation is higher than this value, signal problems and resort to subspaces
		, dumpfileBasename=NULL		##<< scalar string: filename to dump stack before stopping. May set to "recover"
		##end<<
	)
){
	# save calling arguments to allow an continuing an interrupted run
	argsF <- as.list(sys.call())[-1]	# do not store the file name and the first two arguments
	argsF <- argsF[ !(names(argsF) %in% c("","nGen","mc")) ]  # remove positional arguments and arguments mc and nGen
	argsFEval <- lapply( argsF, eval.parent )		# remember values instead of language objects, which might not be there on a repeated call
	#
	#-- fill in default argument values
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
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
		plot( mc0$parms[bo,"a",iPop], mc0$parms[bo,"b",iPop], col=rainbow(100)[twRescale(mc0$resLogDen[bo,"parmsSparce",iPop],c(10,100))] )
	}
	# sum nObs within density
	iDens <- seq_along(mc$dInfos)
	#iDen=1
	iCompsNonFixDen <- lapply( iDens, function(iDen){ 
			irc <-  mc$dInfos[[iDen]]$resCompPos
			irc <- irc[ !(irc %in% iFixTemp) ]
		})
	nObsDen <- sapply( iDens, function(iDen){ sum( nObs[iCompsNonFixDen[[iDen]] ]) })
	TCurr <- getCurrentTemp(mc)
	iBatch=1
	nBatch <- ceiling( nGen/ctrlBatch$nGenBatch )
	for(iBatch in (1:nBatch)){
		if((0 < length(restartFilename)) && is.character(restartFilename) && restartFilename!=""){
			resRestart.twDEMCSA = res #avoid variable confusion on load by giving a longer name
			resRestart.twDEMCSA$iBatch <- iBatch	# also store updated calculation of burnin time
			resRestart.twDEMCSA$args <- argsFEval	# also store updated calculation of burnin time
			save(resRestart.twDEMCSA, file=restartFilename)
			cat(paste("Saved resRestart.twDEMCSA to ",restartFilename,"\n",sep=""))
		}
		resEnd <- thin(res, start=getNGen(res)%/%2 )	# neglect the first half part
		#----- check for high autocorrelation wihin spaces
		mcEnd <- concatPops(stackChainsPop(resEnd))
		#plot( as.mcmc.list(mcEnd), smooth=FALSE )
		specVarRatio <- apply( mcEnd$parms, 3, function(aSamplePurePop){
				spec <- spectrum0.ar(aSamplePurePop)$spec
				varSample <- apply(aSamplePurePop, 2, var)
				spec/varSample
			})	
		if( any(specVarRatio > ctrlConvergence$critSpecVarRatio) ){
			stop("twDEMCSACont: too much autocorrelation. Try using divideTwDEMC")
		} 
		ssc <- stackChainsPop(resEnd)	# combine all chains of one population
		mcl <- as.mcmc.list(ssc)
		#plot( as.mcmc.list(mcl), smooth=FALSE )
		#plot( res$pops[[1]]$parms[,"b",1] ~ res$pops[[1]]$parms[,"a",1], col=rainbow(255)[round(twRescale(-res$pops[[1]]$resLogDen[,"logDen1",1],c(10,200)))] )
		gelmanDiagRes <- try( {tmp<-gelman.diag(mcl); if(length(tmp$mpsrf)) tmp$mpsrf else tmp$psrf[1]} )	# cholesky decomposition may throw errors
		#gelmanDiagRes <- try( gelman.diag(mcl)$mpsrf )	# cholesky decomposition may throw errors
		TEnd <- if(  !inherits(gelmanDiagRes,"try-error") && gelmanDiagRes <= ctrlConvergence$gelmanCrit ){
				logDenT <- calcTemperatedLogDen(resEnd, TCurr)
				#mtrace(getBestModelIndex)
				iBest <- getBestModelIndex( logDenT, resEnd$dInfos )
				maxLogDenT <- logDenT[iBest, ]
				#thetaBest <- stackChains(concatPops(resEnd)$parms)[iBest, ]
				#calcTemperatedLogDen( do.call( dInfos[[1]]$fLogDen, c(list(thetaBest), dInfos[[1]]$argsFLogDen) ), TCurr )
				#sum logDenT within density
				#logDenTDen <- lapply( iDens, function(iDen){
				#			rcpos <- iCompsNonFixDen[[iDen]]
				#			if( length(rcpos)==1) logDenT[,rcpos ,drop=TRUE] else						
				#				rowSums( logDenT[ ,rcpos ,drop=FALSE])
				#		})
				maxLogDenTDen <- sapply(iDens, function(iDen){ 
						sum(maxLogDenT[iCompsNonFixDen[[iDen]] ])
					})
				#TEnd0 <- pmax(1, -2*maxLogDenT/nObs )
				TEndDen <- pmax(1, -2*maxLogDenTDen/nObsDen )
				TEnd <- TCurr
				for( iDen in iDens){
					TEnd[ resEnd$dInfos[[iDen]]$resCompPos ] <- TEndDen[iDen]
				} 
				TEnd[iFixTemp] <- ctrlT$TFix[iFixTemp]
				TEnd[iMaxTemp] <- pmin(TEnd[iMaxTemp], ctrlT$TMax[iMaxTemp])	# do not increase T above TMax
				relTChange <- abs(TEnd - TCurr)/TEnd
				#if( (max(relTChange) <= maxRelTChange) ) recover()
				#trace(isLogDenDrift, recover )
				if( (max(relTChange) <= ctrlConvergence$maxRelTChangeCrit) && !isLogDenDrift(logDenT, resEnd$dInfos, maxDrift=ctrlConvergence$maxLogDenDriftCrit) ){
					res <- resEnd
					print(paste("twDEMCSA: Maximum Temperture change only ",signif(max(relTChange)*100,2),"% and no drift in logDensity. Finishing early.",sep=""))
					break
				}
				# slower TDecrease to avoid Temperatues that are not supported by other datastreams
				TEnd <- TCurr - ctrlT$TDecProp*(TCurr-TEnd)
			}else{
				TEnd <-TCurr
			}
		#TMin <- pmin(TMin, TEnd)
		print(paste("gelmanDiag=",signif(gelmanDiagRes,2)," T=",paste(signif(TCurr,2),collapse=","), sep="") )
		res1 <- res <- twDEMCBlock( resEnd
			, nGen=nGen
			, debugSequential=debugSequential
			, m0=m0
			, controlTwDEMC=controlTwDEMC
			, TEnd=TEnd
			# replace lines below later on by ...
			#, blocks = blocks
			,...
		)
		TCurr <- getCurrentTemp(res)
		print(paste("finished ",iBatch*ctrlBatch$nGenBatch," out of ",nGen," gens.   ", date(), sep="") )
	}
	res$TGlobal <- max(getCurrentTemp(res)[iCompsNonFixDen])
	res$args <- argsFEval
	res
}

.tmp.ggplotResults <- function(){
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
			pred <-  twTwoDenEx1$fModel(ssThin[i,], xSparce=twTwoDenEx1$xSparce, xRich=twTwoDenEx1$xRich, thresholdCovar=thresholdCovar) 
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

.tmp.AllDensities <- function(){ # does not work
	# same as example but with each parameter updated against both densities
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for updating each in turn against both densities
	dInfos=list(
		dSparce=list(fLogDen=denSparcePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		bSparce=list(dInfoPos="dSparce", compPos=c("a","b") )
		,bRich=list(dInfoPos="dRich", compPos=c("a","b") )
	)
	
	do.call( dInfos$dSparce$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparce$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	
	#str(twTwoDenEx1)
	nObs <- c( parmsSparce=1, obsSparce=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	#trace(twDEMCSA, recover )
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs, debugSequential=TRUE, ctrlBatch=list(nGenBatch=256) )
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
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"obsSparce"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"logDen1"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rgb(
			twRescale(logDenT[,"obsSparce"]),0, twRescale(logDenT[,"logDen1"]) ))
	apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )
	
	# density of parameters
	plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")
	
	# predictive posterior (best model only)
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
	
}

.tmp.AllSparce <- function(){ # 
	# same as example but with b also updated against sparce, a only against sparce 
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for updating each in turn against both densities
	dInfos=list(
		dSparce=list(fLogDen=denSparcePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		bSparce=list(dInfoPos="dSparce", compPos=c("a","b") )
		,bRich=list(dInfoPos="dRich", compPos=c("a") )
	)
	
	do.call( dInfos$dSparce$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparce$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	#str(twTwoDenEx1)
	nObs <- c( parmsSparce=1, obsSparce=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	#untrace(twDEMCSA )
	#trace(twDEMCSA, recover )
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs, nBatch=5, debugSequential=TRUE, TFix=c(1,NA,NA) )
	resPops <- res <- twDEMCSACont( res, 256, nObs=nObs, debugSequential=TRUE
		, ctrlT=list( TFix=c(parmsSparce=1) ) 
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
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"obsSparce"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"logDen1"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rgb(
			twRescale(logDenT[,"obsSparce"]),0, twRescale(logDenT[,"logDen1"]) ))
	apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )
	
	# density of parameters
	plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")
	
	# predictive posterior (best model only)
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
	
}


.tmp.2DenSwitched <- function(){ # does not work
	# same as example but with each parameter b updated against sparce
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for updating each in turn against both densities
	dInfos=list(
		dSparce=list(fLogDen=denSparcePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		bSparce=list(dInfoPos="dSparce", compPos=c("b") )
		,bRich=list(dInfoPos="dRich", compPos=c("a") )
	)
	
	do.call( dInfos$dSparce$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparce$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	
	#str(twTwoDenEx1)
	nObs <- c( parmsSparce=1, obsSparce=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
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
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"obsSparce"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"logDen1"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rgb(
			twRescale(logDenT[,"obsSparce"]),0, twRescale(logDenT[,"logDen1"]) ))
	apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )
	
	# density of parameters
	plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")
	
	# predictive posterior (best model only)
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
	
}

.tmp.oneDensity_andCluster <- function(){
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
	nObs <- c( parmsSparce=1, y1=length(twTwoDenEx1$obs$y1), y2=length(twTwoDenEx1$obs$y2) )
	
	#trace(twDEMCSA, recover)
	argsTwDEMCSA <- list( thetaPrior=thetaPrior, covarTheta=covarTheta, dInfos=dInfos, nObs=nObs
		,TFix = c(1,NA,NA)
		,nGen=256
		#,TMax = c(NA,1.2,NA)	# do not increase sparce observations again too much
		, nBatch=3 
		, debugSequential=TRUE
	)
	resPops <- res <- do.call( twDEMCSA, argsTwDEMCSA )
	
	.tmp.continueRun <- function(){
		res$pops[[2]]$spaceInd <- 4		# test spaces different from 1:nSpace
		res2 <- twDEMCSA(res)
	}
	.tmp.byCluster <- function(){
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
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
}


.tmp.den2dCor <- function(){
	# fitting the den2dCor model
	data(den2dCorEx)
	#str3(den2dCorEx)
	
	thetaPrior <- den2dCorEx$thetaPrior
	covarTheta <- diag(den2dCorEx$covarTheta)
	#solve(den2dCorEx$covarTheta)
	invCovarTheta <- 1/covarTheta		# given as independent variances for faster calculation
	
	# for only one logDensity - Temperature of the strongest component goes to zero and other increase according to mismatch
	 dInfos=list( list(fLogDen=den2dCorEx$den2dCor) )
	
	do.call( dInfos[[1]]$fLogDen, c(list(theta=thetaPrior), dInfos[[1]]$argsFLogDen) )
	#str(den2CorEx)
	.nObs <- 1
	
	#trace(twDEMCSA, recover)
	argsTwDEMCSA <- list( thetaPrior=thetaPrior, covarTheta=covarTheta, dInfos=dInfos
		, nObs=.nObs
		,nGen=2*512
		, ctrlBatch=list( 
			nGenBatch=512 
			, useSubspaceAdaptation=TRUE
		)
		, debugSequential=TRUE
		, controlTwDEMC=list(DRgamma=0.1)		# DR step of 1/10 of the proposed lenght
	)
	# if nGen=256, too few samples per space - GelmanDiag does not converge
	res <- res0 <- do.call( twDEMCSA, argsTwDEMCSA )
	res <- res1 <- twDEMCSA( res0, nGen=8*512) 
	res <- res2 <- twDEMCSA( res, nGen=8*512) 
	
	(TCurr <- getCurrentTemp(res))
	#trace(concatPops.twDEMCPops, recover)	#untrace(concatPops.twDEMCPops)
    getSpacesPop(res)
	mc <- concatPops(stackPopsInSpace(res))
	#windows(record=TRUE)
	plot( as.mcmc.list(mc) , smooth=FALSE )
	matplot( mc$temp, type="l" )
	iSpace=1; plot( mc$parms[,"a",iSpace], mc$parms[,"b",iSpace], xlim=c(-0.8,2), ylim=c(-20,20), col=(heat.colors(100))[twRescale(log(-mc$resLogDen[,1,iSpace]),c(10,100))])
	iSpace=2; plot( mc$parms[,"a",iSpace], mc$parms[,"b",iSpace], xlim=c(-0.8,2), ylim=c(-20,20), col=(heat.colors(100))[twRescale(log(-mc$resLogDen[,1,iSpace]),c(10,100))])
	.checkProblemsSpectralPop(res$pops[[2]])
	
	res <- .tmp.byCluster <- function(){
		runClusterParms <- list(
			fSetupCluster = function(){library(twDEMC)}
			,fRun = twDEMCSA
			,argsFRun = within( argsTwDEMCSA,{
					debugSequential<-FALSE 
					#nGen <- 512*2
					nGen <- 512*12
					ctrlConvergence <- list( dumpfileBasename="parms/convergenceProblems" )
				}) 
		)
		# res <- do.call( twDEMCSA, runClusterParms$argsFRun )
		asomRelPath="../../projects/asom"
		save(runClusterParms, file=file.path(asomRelPath,"parms","saDen2dCor.RData"))
		# from asom directory run: 
		#   R CMD BATCH --vanilla '--args paramFile="parms/saDen2dCor.RData"' runCluster.R Rout/runCluster_0.rout 
		# or from libra home directory run 
		#   ./bsubr_i.sh runCluster.R iproc=0 nprocSinge=1 'paramFile="parms/saDen2dCor.RData"'
		# or
		#    bsub -q SLES -n 8 ./bsubr_i.sh runCluster.R  nprocSingle=8 'paramFile="parms/saDen2dCor.RData"'
		load(file.path(asomRelPath,"parms/res_saDen2dCor_0.RData"))
		res <- runClusterRes$res
		#
		res <- loadAssign(file.path(asomRelPath,"parms/restart_saDen2dCor_0.RData"))
		runClusterParms$argsFRun$thetaPrior=res
		save(runClusterParms, file=file.path(asomRelPath,"parms","saDen2dCor.RData"))
		#
		remoteDumpfileBasename = "parms/convergenceProblems"	# search near line 496
		tmp <- loadAssign(file.path(asomRelPath,paste(remoteDumpfileBasename,".rda",sep="")))
		debugger(tmp)
		#	
	 	remoteDumpfileBasename = "parms/res_saDen2dCor_0"	# search near line 496
	 	load(file.path(asomRelPath,paste(remoteDumpfileBasename,".rda",sep="")))
	 	debugger(get(remoteDumpfileBasename))
		
	}
	
}

getBestModelIndex <- function(
	### select the best model based on (temperated) logDensity components
	logDenT		##<< numeric matrix (nStep x nResComp): logDensity (highest are best)
	, dInfos 	##<< list of lists with entry resCompPos (integer vector) specifying the position of result components for each density 
){
	iBest <- if( length(dInfos) > 1){
		# with several densities, each parameter vector is ranked differently
		# select the case where the maximum of all the ranks across densities is lowest
		logDenTDens <- do.call( cbind, 	lapply( seq_along(dInfos), function(iDen){
					dInfo <- dInfos[[iDen]]
					logDenInfoT <- rowSums(logDenT[,dInfo$resCompPos ,drop=FALSE])	
				}))
		rankDen <- apply( -logDenTDens,2, rank) # starting with the highest density
		rankDenMax <- apply( rankDen, 1, max )   
		iBest <- which.min(rankDenMax)
	}else{
		iBest <- which.max( rowSums(logDenT) )
	}
	### the index within logDenT with the best rank
	iBest
}
attr(getBestModelIndex,"ex") <- function(){
	logDenT <- cbind( -sample(5)/2, -sample(5), -sample(5) )
	#dInfos <- list( d1=list(resCompPos=1:2), d2=list(resCompPos=3) )
	dInfos <- list( d1=list(resCompPos=2), d2=list(resCompPos=3) )
	getBestModelIndex(logDenT, dInfos)
	-logDenT
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
	#iDen <- 1
	boL <- sapply( seq_along(dInfos), function(iDen){
			dInfo <- dInfos[[iDen]]
			rs1 <- rowSums(logDenT[1:nr4,  dInfo$resCompPos ,drop=FALSE])
			rs4 <- rowSums(logDenT[(nr+1-nr4):nr, dInfo$resCompPos ,drop=FALSE])
			resTTest <- t.test(rs1,rs4,"less")
			bo <- (diff(resTTest$estimate) >= maxDrift ) && (resTTest$p.value <= alpha) 
		})
	any( boL )
	### TRUE if any of the logDensities are a significantly greater in the fourth quantile compared to the first quantile of the samples
}
attr(isLogDenDrift,"ex") <- function(){
	data(twdemcEx1)
	logDenT <- calcTemperatedLogDen( twdemcEx1, getCurrentTemp(twdemcEx1))
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
.tmp.testRestart <- function(){
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
	else if( nResComp != length( x) ){
		iFix <- sapply( names(x), match, resCompNames )
		if( any(is.na(iFix))) stop("completeResCompVec: x must be of the same length as resCompNames or all its names must correspond names of resCompNames.")
		xAll <- structure( rep( NA_real_, nResComp), names=resCompNames )
		xAll[iFix] <- x
	}
	xAll
	### vector with names resCompNames, with corresponding values of x set, others NA
}







