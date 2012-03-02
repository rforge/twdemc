# see also subspace.R
divideTwDEMCSACont <- function(
	### run twDEMCBlock on subspaces with decreasing Temperature
	mc					##<< former run of twDEMCBlockInt
	, nGen=512			##<< the number of generations for twDEMCBlock
	, nObs							##<< integer vector (nResComp) specifying the number of observations for each result component
	, iPopsDoSplit=integer(0)		##<< populations given here will definitely be splitted (no matter of other criteria)
	#
	, ...							##<< arguments passed to \code{\link{twDEMCBlock}}
	, controlTwDEMC = list()		##<< list argument to \code{\link{twDEMCBlock}} containing entry thin
	, doRecordProposals=FALSE		##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y
	, m0 = calcM0twDEMC(getNParms(mc),getNChainsPop(mc))	##<< minimum number of samples in step for extending runs
	, debugSequential=FALSE			##<< set to TRUE to avoid parallel execution, good for debugging
	, restartFilename=NULL			##<< filename to write intermediate results to
	#
	, ctrlBatch = list(				##<< list of arguments controlling batch executions
		##describe<< 
		#nGenBatch=max(64*controlTwDEMC$thin,m0*(1/ctrlSubspaces$minPSub)*4)		##<< number of generations for one call to twDEMCStep
		nGenBatch=m0*(1/2/ctrlSubspaces$minPSub)*controlTwDEMC$thin*2		##<< number of generations for one call to twDEMCStep
		##<< default: set in a way that on average each population (assuming half are significant) is appended by 2*m0 samples
		, nSampleMin=32				##<< minimum number of samples in each population within batch so that calculation of average density is stable
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
		maxSubPercChangeCrit=1.6	##<< if all subPercChange of all sub-populations are below this value in two last batches, may assume convergence and skip further batches
		, maxRelTChangeCrit=0.025 	##<< if Temperature of the components changes less than specified value, the algorithm can finish
		, maxLogDenDriftCrit=0.3	##<< if difference between mean logDensity of first and fourth quartile of the sample is less than this value, we do not need further batches because of drift in logDensity
		, gelmanCrit=1.4			##<< do not change Temperature, if variance between chains is too high, i.e. Gelman Diag is above this value
		, critSpecVarRatio=10		##<< if proprotion of spectral Density to Variation is higher than this value, split the space
		, pCheckSkipPart=0.5		##<< when checking each population for convergence problems, skip this proportion (kind of burnin) 
		, minNSampleCheck=25		##<< if a population has less samples within a chains, skip the check because it is too uncertain
		, dumpfileBasename=NULL		##<< scalar string: filename to dump stack before stopping. May set to "recover"
	##end<<
	)
	, ctrlSubspaces = list(			##<< list of arguments controlling splitting and merging of subspaces
		##describe<< 
		minPSub = 0.1				##<< minimal proportion of a sub-population
		, maxNSample=max(256,2*m0*getNChainsPop(mc))	##<< if given a value, then looking for subspaces is done on a subsample of given length for efficiency (see \code{\link{getSubSpaces}})
		, argsFSplit=list()			##<< further arguments to \code{\link{findSplit}}
	##end<<
	)
){
	#-- store call arguments 
	# trace(divideTwDEMCSACont, at=3, recover) #untrace(divideTwDEMCSACont)
	# to come past this point as it does not work from within browser	argsF <- as.list(sys.call())[-1]
	argsF <- as.list(sys.call())[-1]	# do not store the file name and the first two arguments
	argsF <- argsF[ !(names(argsF) %in% c("","nGen","mc","iPopsDoSplit")) ]  # remove positional arguments and arguments mc and nGen
	argsFEval <- lapply( argsF, eval.parent )		# remember values instead of language objects, which might not be there on a repeated call
	#
	#-- fill in default argument values
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
	thin <- controlTwDEMC$thin
	frm <- formals()
	ctrlConvergence <- if( hasArg(ctrlConvergence) ) twMergeLists( eval(frm[["ctrlConvergence"]]), ctrlConvergence ) else ctrlConvergence
	ctrlSubspaces <- if( hasArg(ctrlSubspaces) ) twMergeLists( eval(frm[["ctrlSubspaces"]]), ctrlSubspaces ) else ctrlSubspaces 
	ctrlBatch <- if( hasArg(ctrlBatch) ) twMergeLists( eval(frm[["ctrlBatch"]]), ctrlBatch ) else ctrlBatch	# after setting ctrlSubspaces because default uses minPSub
	ctrlT <- if( hasArg(ctrlT) ) twMergeLists( eval(frm[["ctrlT"]]), ctrlT ) else ctrlT
	if( debugSequential ){
		if( 0==length(ctrlConvergence$dumpfileBasename) ) ctrlConvergence$dumpfileBasename <- "recover"
		if( 0==length(ctrlSubspaces$argsFSplit$debugSequential) ) ctrlSubspaces$argsFSplit$debugSequential <- TRUE 
	}
	#
	nResComp <- ncol(mc$pops[[1]]$resLogDen)
	ctrlT$TFix <- completeResCompVec( ctrlT$TFix, colnames(mc$pops[[1]]$resLogDen) )
	iFixTemp <- which( is.finite(ctrlT$TFix) )
	ctrlT$TMax <- completeResCompVec( ctrlT$TMax, colnames(mc$pops[[1]]$resLogDen) )
	iMaxTemp <- which( is.finite(ctrlT$TMax) )	
	#
	ctrlBatch$nGenBatch = (ctrlBatch$nGenBatch%/%thin)*thin		# make nGen multiple of thin
	#boPopsDidConverge <- FALSE
	iGen = 0
	nChainPop <- getNChainsPop(mc)
	#
	iDens <- seq_along(mc$dInfos)
	#iDen=1
	iNonFixTempDens <- lapply( iDens, function(iDen){ 
			irc <-  mc$dInfos[[iDen]]$resCompPos
			irc <- irc[ !(irc %in% iFixTemp) ]
		})
	nObsDen <- sapply( iDens, function(iDen){ sum( nObs[iNonFixTempDens[[iDen]] ]) })
	.spacePop <- getSpacesPop(mc)	# integer vector (nPop): specifying the space replicated that the population belongs to
	spaceInds <- unique(.spacePop)
	.nSamplesPop <- getNSamples(mc)
	.nSamplesSpace <- structure( sapply(spaceInds, function(spaceInd){ sum(.nSamplesPop[.spacePop==spaceInd])}), names=spaceInds)
	#
	m0=calcM0twDEMC(getNParms(mc),getNChainsPop(mc))		##<< minimum number of samples in step for extending runs
	nRowsMin <- max( 2*(m0+1), ceiling(ctrlBatch$nSampleMin/nChainPop) )		# minimum number of rows in population, so that enough samples to split into 2 subs
	nRowsMinMax <- max(nRowsMin, 1/ctrlSubspaces$minPSub*(m0+1))	# must be at least nRowsMin but does not need to be more than samples to split into 1/pSubMin subs
	if( ctrlBatch$nGenBatch%/%thin < nRowsMin) 
		stop("divideTwDEMCSACont: not enough sample per batch. Increase nGenBatch.")
	#
	#-- merge those spaces with less than m0 samples
	pSubs0 <- .nSamplesPop/.nSamplesSpace[ as.character(.spacePop) ]
	resMergePops <- twMergePops(mcApp=mc, m0=m0, minPSub=0, pSubs=pSubs0, nChainPop=nChainPop, spaceInds=spaceInds) # here merge only for m0
	mcApp <- resMergePops$mcApp
	pSubs <- resMergePops$pSubs
	#
	# the following updated variables are required in each cycle
	isMergeAbandonedPops <- FALSE
	mcNewMinN <- subsetF.twDEMCPops(mcApp, fKeep = function(pop){ ### last nRowsMinMax of each chain in each population
			nR <- nrow(pop$parms)
			ret <- max(1, nR+1-nRowsMinMax):nR
		})
	#mcApp <- mcApp		# concatenation of thinned past and new samples
	#pSubs <- pSubs
	subPercChange=rep(ctrlConvergence$maxSubPercChangeCrit, length(pSubs))	# at this value, convergence check is false but pops are splitted (< and >=)	specVarRatioPop <- .checkConvergenceProblems(mcApp, ctrlConvergence$critSpecVarRatio=Inf, pCheckSkipPart=ctrlConvergence$pCheckSkipPart)
	specVarRatioPop <- .checkConvergenceProblems(mcApp, pSubs=pSubs, minPSub=0, nChainPop=nChainPop, pCheckSkipPart=ctrlConvergence$pCheckSkipPart)
	#iPopsDoSplit <- union(iPopsDoSplit, which(specVarRatioPop > ctrlConvergence$critSpecVarRatio))
	#
	while( iGen < nGen ){
		iBatch <- iGen/ctrlBatch$nGenBatch+1
		#if( iBatch > 1){ print("divideTwDEMCSACont: later batch:"); recover() }
		mcApp0 <- mcApp; mcNewMinN0 <- mcNewMinN; pSubs0 <- pSubs	# remember the initial populations
		#mcApp<-mcApp0; mcNewMinN <- mcNewMinN0; pSubs<-pSubs0
		.saveRestartFile( restartFilename, mcApp=mcApp, args=argsFEval )
		#
		#-- can Temperature be decreased?
		TCurr <- getCurrentTemp(mcApp) 
		mcSpace <- stackPopsInSpace(mcApp, mergeMethod="random")	
		mcSpaceEnd <- thin(mcSpace, start=getNGen(mcSpace) %/% 2)
		mcl <- as.mcmc.list(stackChainsPop( mcSpaceEnd ))
		#print("divideTwDEMCSACont: loop start"); recover()
		#plot( mcl, smooth=FALSE )
		#plot( res$pops[[1]]$parms[,"b",1] ~ res$pops[[1]]$parms[,"a",1], col=rainbow(255)[round(twRescale(-res$pops[[1]]$resLogDen[,"logDen1",1],c(10,200)))] )
		gelmanDiagRes <- try( {tmp<-gelman.diag(mcl); if(length(tmp$mpsrf)) tmp$mpsrf else tmp$psrf[1]} )	# cholesky decomposition may throw errors
		TEnd <- .calcTEnd(gelmanDiagRes=gelmanDiagRes, resEnd=mcSpaceEnd, TCurr=TCurr, nObsDen=nObsDen, TFix=ctrlT$TFix, iFixTemp=iFixTemp, iNonFixTempDens=iNonFixTempDens, TMax=ctrlT$TMax, iMaxTemp=iMaxTemp
			, gelmanCrit=ctrlConvergence$gelmanCrit, TDecProp=ctrlT$TDecProp)
		relTChange <- abs(TEnd - TCurr)/TEnd
		#
		#-- check for convergence maybe break
		#print("divideTwDEMCSACont: before checking converge"); recover()
		logDenT <- calcTemperatedLogDen(mcSpaceEnd, TCurr)
		#print("divideTwDEMCSACont: before checking convergence"); recover()
		if( 
			gelmanDiagRes <= ctrlConvergence$gelmanCrit &&		
			(maxRelTChange <- max(relTChange)) <= ctrlConvergence$maxRelTChangeCrit && 
			(maxspecVarRatioPop <- max(specVarRatioPop)) < ctrlConvergence$critSpecVarRatio &&		
			# important not <= because of initialization (do not break before first batch where subPercChange is unknown)
			(maxSubPerChange <- max(subPercChange[pSubs >= ctrlSubspaces$minPSub])) < ctrlConvergence$maxSubPercChangeCrit &&	
			!isLogDenDrift(logDenT, mcSpaceEnd$dInfos, maxDrift=ctrlConvergence$maxLogDenDriftCrit)  
			){
			cat(paste("divideTwDEMCSACont: Finishing early."
					,"\n  gelman.diag =",signif(gelmanDiagRes,2)
					,"\n  max spectral density to variance ratio =",signif(maxspecVarRatioPop,2)
					,"\n  max change of subspaces proportions=",signif(maxSubPerChange,2)
					,"\n  max Temperture change =",signif(maxRelTChange*100,2),"%"
					,"\n  and no drift in logDensity"
					,"\n",sep=""))
			break
		}
		#
		#-- split pops that changed poportion a lot or pops with high spectralDensity ratio into smaller pops
		# important > instead of >= maxSubPercChangeCrit; because of initialization: subPercChange is unknown
		iPopsSplit <- union( iPopsDoSplit, which( (pSubs/2 > ctrlSubspaces$minPSub) & 
					((subPercChange > ctrlConvergence$maxSubPercChangeCrit) | (specVarRatioPop > 0.8*ctrlConvergence$critSpecVarRatio)) 
			))
		#.getParBoundsPops(mcApp$pops)
		#.getParBoundsPops(mcAppRecent$pops)
		#trace(.splitPops, recover )
		resSplitPops <- .splitPops(mcNewMinN=mcNewMinN, mcApp=mcApp, pSubs
			#, subPercChange=subPercChange, specVarRatioPop=specVarRatioPop, critSpecVarRatio=ctrlConvergence$critSpecVarRatio
			#,iPopsDoSplit=iPopsDoSplit, subPercChangeCrit=subPercChangeCrit
			,iPopsSplit=iPopsSplit
			#,minPSub=ctrlSubspaces$minPSub, maxNSample=ctrlSubspaces$maxNSample
			, m0=m0, nChainPop=nChainPop, spaceInds=spaceInds
			, ctrlSubspaces=ctrlSubspaces
		#,argsFSplit=ctrlSubspaces$argsFSplit	##<< further arguments to \code{\link{findSplit}}
		)
		#if( any(getNSamples(resSplitPops$mcApp)<m0)) stop("divideTwDEMCSteps (afterSplit): encountered too few samples in population.")
		mcApp <- mcApp1 <- resSplitPops$mcApp
		mcNewMinN <- mcNewMinN1 <- resSplitPops$mcNewMinN
		pSubs <- pSubs1 <- resSplitPops$pSubs
		iPopsDoSplit <- integer(0)
		#
		#-- merge subspaces that contain less samples than minPSub/2
		#print("divideTwDEMCSACont: before merge"); recover()
		#trace(twMergePops, recover)	#untrace(twMergePops)
		if( isMergeAbandonedPops ){
			resMergePops <- twMergePops(mcApp=mcApp, pSubs=pSubs, minPSub=ctrlSubspaces$minPSub, nChainPop=nChainPop, m0=m0, spaceInds=spaceInds)
			mcApp <- mcApp3 <- resMergePops$mcApp
			pSubs <- pSubs3 <- resMergePops$pSubs
			mcNewMinN <- mcNewMinN3 <- resMergePops$mcNewMinN	# TODO: merge also the mcNewMinN
			subPercChange <- subPercChange[resMergePops$iPopsBefore]
			specVarRatioPop <- specVarRatioPop[resMergePops$iPopsBefore]
		}
		#
		#-- sample the next batch
		print(paste("gelmanDiag=",signif(gelmanDiagRes,2)," T=",paste(signif(TCurr,2),collapse=","), sep="") )
		if( any( (.nsMinN <-getNSamples(mcNewMinN)) < m0) ){
			if( any(.nsMinN < m0/2) )
				warning(paste("divideTwDEMCSACont: Very few samples after splitting mcNew: m0=",m0," nSamples=",paste(.nsMinN[.nsMinN<m0],collapse=","),sep="")); 
			print(paste("divideTwDEMCSACont: May indicate problems: Fewer samples than m0 after splitting mcNew: m0=",m0," nSamples=",paste(.nsMinN[.nsMinN<m0],collapse=","),sep="")); 
			# may occur because minPSub is done on a subSample in .sp
		}
		nGenStep <- min(ctrlBatch$nGenBatch, nGen-iGen)
		#mtrace(divideTwDEMCStep)
		#resStep <- divideTwDEMCStep(mcApp, nGen=nGenStep, doRecordProposals=doRecordProposals, m0=m0, TEnd=TEnd, ... )
		# initialized with mcNewMinN, because some pops of mcApp might have too few samples
		resStep <- divideTwDEMCStep(mcNewMinN, qPop=pSubs, nGen=nGenStep, doRecordProposals=doRecordProposals, m0=m0, nRowsMin=nRowsMin, nRowsMinMax=nRowsMinMax, minPSub=ctrlSubspaces$minPSub, TEnd=TEnd, ... )
		mcNew <- mcNew2 <- resStep$resTwDEMC			# mcNew: only the new samples
		print(paste("nGen(mcNew)=",paste(getNGen(mcNew),collapse=",")))		
		mcNewMinN <- mcNew2 <- resStep$resTwDEMCMinN	# mcNew: last 2*m0 new rows 
		pSubs <- pSubs2 <- pSubsNew <- resStep$pSubs	
		subPercChange <- resStep$subPercChange
		.spacesPop <- getSpacesPop(mcNew)
		if( !isTRUE( all.equal(.pSpaces <- as.numeric(tapply(pSubs, .spacesPop, sum)), rep(1,length(unique((.spacesPop)))) )) )
			stop("divideTwDEMCSteps (1): pSubs within subspaces do not sum to 1.")
		#
		#-- check convergence problems of the new samples
		#specVarRatioPop <- .checkConvergenceProblems(mcNew, critSpecVarRatio=Inf, pCheckSkipPart=pCheckSkipPart)
		#specVarRatioPop <- .checkConvergenceProblems(mcNew, pSubs=pSubs, minPSub=ctrlSubspaces$minPSub, nChainPop=nChainPop, critSpecVarRatio=ctrlConvergence$critSpecVarRatio, dumpfileBasename=ctrlConvergence$dumpfileBasename, pCheckSkipPart=ctrlConvergence$pCheckSkipPart, minNSampleCheck=ctrlConvergence$minNSampleCheck)
		specVarRatioPop <- .checkConvergenceProblems(mcNewMinN, pSubs=pSubs, minPSub=ctrlSubspaces$minPSub, nChainPop=nChainPop, critSpecVarRatio=ctrlConvergence$critSpecVarRatio, dumpfileBasename=ctrlConvergence$dumpfileBasename, pCheckSkipPart=ctrlConvergence$pCheckSkipPart, minNSampleCheck=ctrlConvergence$minNSampleCheck)
		#
		#-- append to thinned past
		#.tmp <- lapply( mcNew$pops[nSamplePop != 0], .checkPop, nGen=12 )		
		#.spacesPop <- getSpacesPop(mcNew)
		#tapply(nSamplePop, .spacesPop, function(nsI){ nsI/sum(nsI)})
		#
		#.spacesPopOrig <- getSpacesPop(mcApp)
		#.pSpacesOrig <- as.numeric(tapply(pSubs, .spacesPopOrig, sum)
		#boPopsDistConverge <- all(resStep$subPercChange <= subPercChangeCrit)
		mcAppPast <- squeeze(mcApp, length.out=ceiling(getNSamples(mcApp)*ctrlBatch$pThinPast) ) # mcApp0: thinned former mcApp
		mcApp <- mcApp2 <- .concatTwDEMCRuns(mcAppPast,mcNew,doRecordProposals=doRecordProposals) # mcApp and mcApp1: thinned past + new samples 
		#.nS1 <- sum(getNSamples(mcApp2))
		#if( nGenStep/thin >= 64){  # only check and split populations, if enough samples in batch
		#iPop=length(mcNew$pops)
		#
		#} # enough new samples for checking/splitting
		#
		iGen <- iGen + nGenStep
		#cat(paste(iGen," out of ",nGen," generations completed. T=",paste({T<-res$temp[nrow(res$temp),];round(T,digits=ifelse(T>20,0,1))},collapse="  "),"     ",date(),"\n",sep=""))
		iPopsLarge <- which( pSubs >= quantile(pSubs,0.4))	# neglect the small populations in subPercChange calculation 
		maxSubPercChange <- max(subPercChange[iPopsLarge])
		cat(paste(iGen," out of ",nGen," gens completed. ,max(subPercChange)=",signif(maxSubPercChange,2),"   nPops=",paste(resStep$nPopsApp,collapse=","),"     ",date(),"\n",sep=""))
	}	# while iGen < nGenBatch
	##value<< twDEMC with additional list entries
	resAdd <- list(
		nGen = iGen						##<< scalar integer: number of completed generations
		,pSubsBatch = pSubs					##<< numeric vector (nPop): proportion of the populations within space during last batch
		,subPercChange = subPercChange		##<< numeric vector(nPop): relative change of proportions during last batch 
		,relTChange =  relTChange			##<< numeric vector (nDen): calculated relative change of calculated new temperature for the next batch 
		,specVarRatioPop = specVarRatioPop	##<< numeric vector (nPop): ratio of spectral density to variance
		,TGlobal = max(getCurrentTemp(mcApp)[iNonFixTempDens]) ##<< Temperature of non-fixed components
		,args=argsFEval						##<< calling arguments to provide restart capability
	##end<<
	)
	mcApp[ names(resAdd) ] <- resAdd
	mcApp
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
	nObs <- 1
	#
	#trace(twDEMCSA, recover)
	argsTwDEMCSA <- list( thetaPrior=thetaPrior, covarTheta=covarTheta, dInfos=dInfos, nObs=nObs
		, nGen=256
		, nBatch=1 
		, debugSequential=TRUE
		, controlTwDEMC=list(DRgamma=0.1)		# DR step of 1/10 of the proposed lenght
	)
	resPops <- res <- res0 <- do.call( twDEMCSA, argsTwDEMCSA )
	res$pops[[2]]$spaceInd <- 4
	#trace(twMergePops, recover)	#untrace(twMergePops)
	#trace(divideTwDEMCSACont, at=3, recover)	#untrace(divideTwDEMCSACont)
	resDiv <- res2Div <- divideTwDEMCSACont( res, nObs=nObs, nGen=1024 
		, debugSequential=TRUE
		, controlTwDEMC=list(DRgamma=0.1)		    # DR step of 1/10 of the proposed length
		, iPopsDoSplit=1:2			# definitely split the populations
		, ctrlSubspaces=list(maxNSample=250)
	)
	#str(resDiv$args)
	resDiv <- res3Div <- do.call( divideTwDEMCSACont, c(list( res2Div, 1*1024), resDiv$args) )
	resDiv <- res3Div <- do.call( divideTwDEMCSACont, c(list( resDiv, 4*1024), twMergeLists( resDiv$args, list( 
					ctrlConvergence=list(gelmanCrit=1.2)
				))))
	#
	resDiv$pSubs
	resDiv$subPercChange
	getNSamples(resDiv)
	res1 <- resDiv
	res2 <- stackChainsPop(tmp <- stackPopsInSpace(resDiv, mergeMethod="random", nInSlice = 1))
	mc <- thin(concatPops(res2), start=256)
	#plot( tail(tmp$pops[[1]]$temp) )
	#plot( res2$pops[[1]]$temp )
	getCurrentTemp(res1)
	plot( as.mcmc.list(res2), smooth=FALSE )
	#.getParBoundsPops(resDiv$resTwDEMC$pops)
	#trace(subset.twDEMC,recover)	#untrace(subset.twDEMC)
	matplot( mc$pAccept[,1,], type="l" )
	matplot( mc$temp, type="l" )
	matplot( mc$resLogDen[,1,], type="l" )
	iSpace=2; plot( mc$parms[,"a",iSpace], mc$parms[,"b",iSpace], xlim=c(-0.8,2), ylim=c(-20,20), col=(heat.colors(100))[twRescale(log(-mc$resLogDen[,1,iSpace]),c(10,100))])
	gelman.diag(as.mcmc.list(res2))
	#trace(.checkProblemsSpectralPop,recover)	#untrace(.checkProblemsSpectralPop)
	.checkProblemsSpectralPop(resDiv$pops[[1]])
	effectiveSize(as.mcmc.list(res2))
	
	# no mixing, stays at a given Temperature
	# detects autocorrelation
	TCurr <- getCurrentTemp(resPops)
	set.seed(0815)
	resd <- divideTwDEMCSteps(resPops
		, nGen=256*5
		, nGenBatch=256
		, dInfos=dInfos
		, debugSequential=TRUE
		, controlTwDEMC=list(DRgamma=0.1)
		, minPSub=0.05
		, TSpec=cbind( T0=TCurr, TEnd=TCurr )
	)
	str3(resd)
	plot( as.mcmc.list(stackPopsInSpace(resd$resTwDEMC)), smooth=FALSE )
	plot( as.mcmc.list(stackChainsPop(stackPopsInSpace(resd$resTwDEMC))), smooth=FALSE )
	resd$subPercChange
}



divideTwDEMCStep <- function(
	### run twDEMCBlock on subspaces
	aTwDEMC					##<< former run of twDEMCBlockInt
	,qPop=numeric(0)		##<< numeric vector (nPop): probability of samples in subspace within entire space 
	,nGen=512				##<< the number of generations for twDEMCBlock
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMCBlock}} containing entry thin
	, m0 = calcM0twDEMC(getNParms(aTwDEMC),getNChainsPop(aTwDEMC))	##<< minimum number of samples in step for extending runs
	, nRowsMin				##<< number of rows in sub to sample at minimum, must contain enough samples for robust estimation of mean log-Density of subspace
	, nRowsMinMax			##<< number of rows in sub to include as maximum mcNewMinN
	, minPSub				##<< minimum proportion of a subsample in splitting´
	, TEnd					##<< numeric vector (nResComp) target temperature 
	, ...					##<< further arguments to \code{\link{twDEMCBlock}}, such as TEnd
){
	if( length(aTwDEMC$dInfos) > 1)
		warning("divideTwDEMCStep: subspace log-Density weighted aggregation only possible with one log-density function.")
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
	thin <- controlTwDEMC$thin
	nGenThin <- (nGen %/% thin)*thin		# when dealing with thinned samples use the number of generations of last recorded thinned generation
	nBlock <- getNBlocks(aTwDEMC)		# number of blocks in aTwDEMC
	nPop <- getNPops(aTwDEMC)			# number of populations
	iPops <- 1:nPop						
	spacesPop <- getSpacesPop(aTwDEMC)
	spaceInds <- unique(spacesPop)
	nPopSpace <- table(spacesPop)		# number of populaions per space replicated
	nSpace <- length(nPopSpace)			# number of space replicates
	iPopsSpace <- lapply(spaceInds, function(iSpace){ which(spacesPop==iSpace)}) # pops in flat version per space
	names(iPopsSpace) <- spaceInds
	nSamplePop <- getNSamples(aTwDEMC)
	#
	#----  calculating initial quantiles and number of generations
	if( 0 == length(qPop) ){
		#iiSpace=nSpace
		qPopSpace <- lapply( 1:nSpace, function(iiSpace){
				nSampleIPop <- nSamplePop[iPopsSpace[[iiSpace]] ]
				nSampleIPop / sum(nSampleIPop)
			} ) 
		qPop <- do.call(c,qPopSpace)
	}
	#nGen0 <- pmax(m0*thin, ceiling(minNSamplesSub*thin/nChainPop), nGen*qPop)		# at minimum m0*thin generations to keep sufficient samples for extending the run
	nGen0Pops <- nGen*qPop
	#nGen0PopsThin <- pmax( m0, nSampleMin, ceiling(nGen0Pops/thin))*thin	# must have at least m0 samples in order to be able to extend the sample
	nGen0PopsThin <- rep( nRowsMin*thin, nPop )
	#
	#----- mark the pops, for which to append samples
	iPopsApp <- which(qPop >= 0.9*minPSub)
	#iiSpace <- nSpace
	for( iiSpace in 1:nSpace){
		iSpace <- spaceInds[iiSpace]
		iPopsAppSp <- iPopsApp[ spacesPop[iPopsApp] == iSpace ]
		iPopsNoAppSp <- setdiff( iPopsSpace[[iiSpace]], iPopsAppSp )
		# also append pops that have a common border with the pops of significant sample proportions
		if( length(iPopsNoAppSp) != 0){
			# extract lower and upper ParBounds values for each parameter
			upb <- do.call( c ,lapply( aTwDEMC$pops[iPopsAppSp], "[[", "upperParBounds"))
			lpb <- do.call( c, lapply( aTwDEMC$pops[iPopsAppSp], "[[", "lowerParBounds"))
			upbV <- lapply( unique(names(upb)), function(pName){ upb[names(upb)==pName] })
			names(upbV) <- unique(names(upb))
			lpbV <- lapply( unique(names(lpb)), function(pName){ lpb[names(lpb)==pName] })
			names(lpbV) <- unique(names(lpb))
			#iPop <- iPopsNoAppSp[1]
			for( iPop in iPopsNoAppSp){
				#  check if a dimension of lowerParBound matches an upperBarBound of pops that need to be appended
				lpb2 <- aTwDEMC$pops[[iPop]]$lowerParBounds
				#pName <- names(lpb2)[1]
				boCommonBorder = FALSE
				for( pName in names(lpb2) ){ 
					if( lpb2[pName] %in% upbV[[pName]]){
						boCommonBorder <- TRUE
						break
					} 
				}#for lpb2
				if( !boCommonBorder ){
					upb2 <- aTwDEMC$pops[[iPop]]$lowerParBounds
					#pName <- names(lpb2)[1]
					for( pName in names(upb2) ){ 
						if( upb2[pName] %in% lpbV[[pName]]){
							boCommonBorder <- TRUE
							break
						} 
					}#for lpb2
				}
				if( boCommonBorder ) iPopsApp <- c(iPopsApp, iPop)
			} # for iPop
		} # if there are abandoned Pops
	} # for iiSpace
	iPopsNoApp <- setdiff( seq_along(qPop), iPopsApp )
	nGen0PopsThin[iPopsNoApp] <- 0
	#
	#----- initial run based on current quantiles
	#trace("twDEMCBlockInt")
	resTwDEMC0 <- twDEMCBlock(aTwDEMC, nGen=nGen0PopsThin, extendRun=FALSE, controlTwDEMC=controlTwDEMC, m0=m0, TEnd=TEnd, ...)
	# keep the samples of the abandoned spaces, but set temperature of the last row to TEnd
	resTwDEMC0$pops[ iPopsNoApp ] <- aTwDEMC$pops[ iPopsNoApp ]	
	for( iPop in iPopsNoApp){
		resTwDEMC0$pops[[iPop]]$temp[nrow(resTwDEMC0$pops[[iPop]]$temp),] <- TEnd
	}
	resTwDEMC <- resTwDEMC0
	#XX Think about temperature
	
	#plot(as.mcmc.list(resTwDEMC),smooth=FALSE)
	#mcp <- concatPops( mcpBlock <- subPops(resTwDEMC0, iPops=nPop))
	#plot( as.mcmc.list(mcp), smooth=FALSE)
	#
	#----- calculating quantiles based on density of the (pSubs) 
	##details<< 
	## A first estimate of the proportion of samples from different subspaces
	## are the initial percentiles qp. 
	## The proportion of the samples from different subspaces is estimated
	## by the proportions of integrated density of the subspaces.
	## These proportions are esimated by average density over all samples multiplied by an
	## estimate proportional to the volume: the initial quantiles.
	#iPop <- 1
	#iPop <- 2
	#print("divideTwDEMCStep: Before calculating mean logDen."); recover()
	popLogMeanDensSubs <- lapply( 1:nSpace, function(iiSpace){
			iPops <- iPopsSpace[[iiSpace]]
			#iPop <- iPops[length(iPops)]
			maxLogDenPops <- sapply(iPops,function(iPop){max( resTwDEMC$pops[[iPop]]$resLogDen, na.rm=TRUE )})
			maxLogDen <- max(maxLogDenPops)
			lw <- sapply(iPops,function(iPop){  # average the unnormalized densities
					#tempLogDen <- calcTemperatedLogDenChains( resTwDEMC$pops[[iPop]]$resLogDen, TEnd)
					#logDen <- sumLogDenComp(tempLogDen, resTwDEMC$dInfos)
					#assume only one density
					tempLogDen <- calcTemperatedLogDen( stackChains(resTwDEMC$pops[[iPop]]$resLogDen), TEnd, maxLogDen=maxLogDen )
					#ss <- rowSums(tempLogDen)
					#ss <- stackChains(resTwDEMC$pops[[iPop]]$logDen)
					ssLogDen <- rowSums(tempLogDen)
					twLogSumExp(ssLogDen, shiftUpperBound=TRUE)-log(length(ssLogDen)) 
				})	
			lw - max(lw, na.rm=TRUE)			#divide by a constant on exp scale to avoid numerical errors
		})	
	logMeanDensSubs <- do.call(c, popLogMeanDensSubs )
	#barplot(logMeanDensSubs)
	#
	# estimate proportion of subspaces in the limiting distribution
	# it will be used to sample from the subspaces 
	# iiSpace=iSpaces[1]
	pPops <- numeric( length(spacesPop))
	for( iiSpace in 1:nSpace){
		p2u <- exp(popLogMeanDensSubs[[iiSpace]])*qPop[iiSpace]		# two weights: quantile and density importance
		pPops[ iPopsSpace[[iiSpace]] ] <- p2 <- p2u/sum(p2u)									# normalize to 1 within population 
	}
	subPercChange <- pPops/qPop	# ratio of estimated proportion in limiting distribution to proportion of initial proportion
	nSamplesSubsReq0 <- nGenThin/thin * pPops
	# because we sampled at least nSampleMin, we can append a factor of required samples
	#
	# because of rounding small differences in sum of sample numbers may occur 
	# for small deviations, adjust the number of samples in the largest sample
	# add argument nGen to pops
	#iiSpace=nSpace
	nSamplesSubsReq <- rep( NA_integer_, length(nSamplesSubsReq0) )
	for( iiSpace in 1:nSpace ){ 
		extFac <- max(1,floor(min( (nGen0PopsThin/thin /nSamplesSubsReq0)[ iPopsSpace[[iiSpace]] ])))
		nSamplesSubsReq[ iPopsSpace[[iiSpace]] ] <- nsReq <- round(extFac*nGenThin/thin * pPops[ iPopsSpace[[iiSpace]] ])	# at least one sample in the subspace
	}
	#dGen <- extFac*nGenThin/thin - sum(nsReq)
	#	if( dGen != 0 ){
	#			if( abs(dGen)>length(iPopsSpace[[iiSpace]]) ) stop("divideTwDEMC: problem in calculating subspace sample numbers")
	#		iSubMax <- iPopsSpace[[iiSpace]][ which.max( nSamplesSubsReq[iPopsSpace[[iiSpace]] ]) ]
	#		nSamplesSubsReq[iSubMax] <- nSamplesSubsReq[iSubMax] + dGen
	#	} 
	#}
	#nSamplesSubsReq <- pmax(1,nSamplesSubsReq)	# keep at least 1 sample per subspace
	#
	#-------  extend the previous runs
	nSamples0 <- getNSamples(resTwDEMC0)
	nSamplesAdd <- pmax(0, nSamplesSubsReq-nSamples0)
	resTwDEMC <- resTwDEMC1 <- if( max(nSamplesAdd) == 0) resTwDEMC0 else{
			# stay at current temperature (we are already at the prescribed end temperature)
			#for( iPop in seq_along(resTwDEMC0$pops) ){
			#	resTwDEMC0$pops[[iPop]]$TEnd <- 
			#		resTwDEMC0$pops[[iPop]]$temp[ nSamples0[iPop] ]
			#}
			#mtrace(twDEMCBlock.twDEMCPops)
			#mtrace(twDEMCBlockInt)
			resTwDEMC <- twDEMCBlock( resTwDEMC0, nGen=nSamplesAdd*thin, controlTwDEMC=controlTwDEMC,  m0=m0, TEnd=TEnd, ...  )
		}
	#all( getNSamples(resTwDEMC) >= nSamplesSubsReq)
	#getNSamples(resTwDEMC1)
	
	#----- do a subsampling
	# because we required minNSamples to take we possibly can add more than required
	resTwDEMC <- resTwDEMC2 <- squeeze.twDEMCPops(resTwDEMC1, length.out=nSamplesSubsReq)
	#getNSamples(resTwDEMC2)
	# set Temperatreu to end Temperature

	#----- also provide unthinned last 2*m0 samples of each population
	# when providing only the thinned sample, some subs may have too few samples to split and restart
	# assumes integer nRowsMinMax > 0
	#print("divideTwDEMCStep: before subsetting last nRowsMin"); recover()
	fKeep <- function(pop){ ### last 2*m0 rows of each chain 
		nR <- nrow(pop$parms)
		ret <- max(1, nR+1-nRowsMinMax):nR
	}
	#trace(subsetF.twDEMCPops, recover) 	#untrace(subsetF.twDEMCPops)
	resTwDEMCMinN <- subsetF.twDEMCPops(resTwDEMC1, fKeep)
	#tail(resTwDEMCMinN$pops[[4]]$temp)

	##value<< 
	## For each population, a list with entries
	res <- list(	##describe<<
		resTwDEMC = resTwDEMC	##<< the runs with new samples of class twDEMCPops
		,pSubs = pPops			##<< the quantiles of samples within given subspace
		,subPercChange = subPercChange	##<< numeric vector: ratio of estimated proportion in limiting distribution to proportion of initial proportion
		,resTwDEMCMinN = resTwDEMCMinN	##<< twDEMCPops with last (nRowsMin <- max( 2*m0, nSampleMin )) rows
		,nPopsApp = table( spacesPop[iPopsApp] )	##<< numeric vectgor (nSpace): number of spaces for which samples are appended
		)
	##end<<
}
attr(divideTwDEMCStep,"ex") <- function(){
	data(den2dCorEx)
	aTwDEMC <- den2dCorEx$mcSubspaces0
	getSpacesPop(aTwDEMC)
	getNSamples(aTwDEMC)

	#mtrace(divideTwDEMCStep)
	res <- divideTwDEMCStep(aTwDEMC, nGen=256, dInfos=list(list(fLogDen=den2dCor))
		,  debugSequential=TRUE
		#,  controlTwDEMC=list(DRGamma=0.05) 
	)
	getNSamples(res$resTwDEMC)
	
	#windows(record=TRUE)
	plot( as.mcmc.list(stackPopsInSpace(aTwDEMC)), smooth=FALSE)
	#mc1 <- stackPopsInSpace(res$resTwDEMC, mergeMethod="stack")
	mc1 <- stackPopsInSpace(res$resTwDEMC, mergeMethod="random")
	plot( as.mcmc.list(mc1), smooth=FALSE) # note how the distribution shifted
	#plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=4)), smooth=FALSE)
	#plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=30)), smooth=FALSE)
	
	ss <- stackChains(concatPops(mc1))			# the new sample
	ss0 <- stackChains(concatPops(squeeze(stackPopsInSpace(aTwDEMC),length.out=getNSamples(mc1)[1]))) # initial sample of same length
	plot( b ~ a, as.data.frame(ss0), xlim=c(-0.5,2), ylim=c(-20,40) )
	plot( b ~ a, as.data.frame(ss), xlim=c(-0.5,2), ylim=c(-20,40) ) # not that more samples are in the region of interest
	
	.tmp.f <- function(){
		#barplot(res[[1]]$subPercChange, names.arg=seq_along(res[[1]]$subPercChange) )
		#str(res[[1]]$lowerParBounds)
		ssImpPops1 <- ssImpPops <- abind( lapply( res, "[[", "sample"), rev.along=0 )
		plot(density( ssImpPops[,"a",1]));lines(density( ssImpPops[,"a",2]),col="green"); lines(density( ss0[,"a"]),col="blue")
		plot( b ~ a, as.data.frame(ss0), xlim=c(-0.5,2), ylim=c(-20,40) ); points(0.8, 0, col="red" )
		#plot( b ~ a, as.data.frame(ssImpPops[,,2]) ); points(xyMax[1], xyMax[2], col="red" )
		plot( b ~ a, as.data.frame(ssImpPops[,,2]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(xyMax[1], xyMax[2], col="red" )
		plot( b ~ a, as.data.frame(ssImpPops[,,1]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(xyMax[1], xyMax[2], col="red" )
		#plot(density( ssImpPops[,"b",1]));lines(density( ssImpPops[,"b",2]),col="green"); lines(density( ss0[,"b"]),col="blue")
		ssImpPops2 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops1[,,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE ), "[[", "sample"), rev.along=0 )
		#mtrace(divideTwDEMC)
		#ssImpPops2 <- ssImpPops <- divideTwDEMC(ssImpPops1[,-1,], nGen=100, fLogDen=den2dCor, attachDetails=TRUE )
		ssImpPops3 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops2[,,], nGen=500, fLogDen=den2dCor ), "[[", "sample"), rev.along=0 )
		ssImpPops4 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops3[,,], nGen=500, fLogDen=den2dCor ), "[[", "sample"), rev.along=0 )
		ssImpPops5 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops4[,,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE ), "[[", "sample"), rev.along=0 )
		ssImpPops6 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops5[,,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE ), "[[", "sample"), rev.along=0 )
		str(ssImpPops)
	}
}


.checkProblemsGelman1 <- function(
	### convergence between several populations
	aTwDEMC
	, criticalGelmanDiag=1.2
){
	
	resCoda <- as.mcmc.list(aTwDEMC)
	#str(resCoda)
	#plot( resCoda, smooth=FALSE )
	tmp2 <- gelman.diag(resCoda)
	list(
		hasProblem = (tmp2$mpsrf > criticalGelmanDiag)
		,resGelmanDiag = tmp2
	)
}

.checkProblemsSpectralPop <- function(
	### calculating the spectral density for each chain of pop
	pop						##<< list: the population to check
	, critSpecVarRatio=25	##<< If the ratio of the spectralDensity to Variance is above this value, stop 
){
	#parmsChain <- adrop(pop$parms[,,1 ,drop=FALSE],3)
	# 
	ss <- stackChains(pop$parms)
	spec <- try( spectrum0.ar(ss)$spec, silent=TRUE )
	if( !inherits(spec,"try-error") ){
		varSample <- apply(ss, 2, var)
		specVarRatioLumped <- ifelse( varSample==0, Inf, spec/varSample)
	}else{
		return( list(
			hasProblems=TRUE
			,specVarRatioLumped = Inf
				))
	}
	if( all(specVarRatioLumped < critSpecVarRatio) ){
		return( list(
			hasProblems=FALSE
			,specVarRatioLumped = specVarRatioLumped
		) )
	}else{
		specVarRatio <- apply( pop$parms, 3, function(parmsChain){
				spec <- spectrum0.ar(parmsChain)$spec
				varSample <- apply(parmsChain, 2, var)
				ifelse( varSample==0, Inf, spec/varSample)
			})
		list(
			hasProblems = TRUE
			,specVarRatioLumped = specVarRatioLumped
			,specVarRatio = specVarRatio 
		)
	}
}

.tmp.f <- function(){}



.saveRestartFile <- function(
	### saving a restart file
	restartFilename	##<< the filename
	,mcApp			##<< the object to save
	,args=list()	##<< list of arguments for the restart
	,msg="Saved resRestart.twDEMCSA to "	##<< printout msg prepended to filename, set to NULL to prevent output
){
	if((0 < length(restartFilename)) && is.character(restartFilename) && restartFilename!=""){
		resRestart.twDEMCSA = mcApp #avoid variable confusion on load by giving a longer name
		resRestart.twDEMCSA$args <- args		# also store updated calculation of burnin time
		save(resRestart.twDEMCSA, file=restartFilename)
		if( 0 != length(msg) ) cat(paste(msg,restartFilename,"\n",sep=""))
	}
}

.calcTEnd <- function(
	### calculating target temperature
	gelmanDiagRes	##<< gelman criterion
	,resEnd			##<< twDEMC results with concatenated spaces 
	,TCurr			##<< numeric vector (nResComp): current Temperature
	,nObsDen		##<< integer vector (nDens): number of observations for each density
	,TFix			##<< numeric vector (nResComp): fixed temperature for several components, other NA
	,iFixTemp=which(is.finite(TFix))	##<< positions of result components with fixed density
	,iNonFixTempDens 	##<< list (iDen): with integer vector entries: positions of result components that have no fixed density for each density
	,TMax			##<< numeric vector (nResComp): maximum temperature for several components, other NA
	,iMaxTemp=which(is.finite(TMax))	##<< positions of result components with maximum Temperature
	, TDecProp=0.9				##<< proportion of Temperature decrease: below one to diminish risk of decreasing Temperature too fast (below what is supported by other data streams)
	, gelmanCrit=1.4			##<< do not decrease Temperature, if variance between chains is too high, i.e. Gelman Diag is above this value
){
	#print(".calcTEnd: start"); recover()
	#
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
	iDens <- seq_along(iNonFixTempDens)
	maxLogDenTDen <- sapply(iDens, function(iDen){ 
			sum(maxLogDenT[iNonFixTempDens[[iDen]] ])
		})
	#TEnd0 <- pmax(1, -2*maxLogDenT/nObs )
	TEndDen <- pmax(1, -2*maxLogDenTDen/nObsDen )
	TEnd <- TCurr
	for( iDen in iDens){
		TEnd[ resEnd$dInfos[[iDen]]$resCompPos ] <- TEndDen[iDen]
	} 
	TEnd[iFixTemp] <- TFix[iFixTemp]
	TEnd[iMaxTemp] <- pmin(TEnd[iMaxTemp], TMax[iMaxTemp])	# do not increase T above TMax
	# slower TDecrease to avoid Temperatues that are not supported by other datastreams
	TEnd <- ifelse( TEnd < TCurr, TCurr - TDecProp*(TCurr-TEnd), TEnd )
	#
	if(  inherits(gelmanDiagRes,"try-error") || gelmanDiagRes > gelmanCrit ){		
		TEnd <-pmax(TEnd, TCurr)
	}
	### numeric vector (nResComp): target temperature for all result components
	if( any(TEnd < 1)) stop("encountered TEnd < 1")
	TEnd
}

.checkConvergenceProblems <- function(
	### checking all population for convergence problems
	mc			##<< result of twDEMC, usually only the last samples without the history
	,pSubs						##<< current percentiles of population
	, minPSub = 0.1				##<< minimal quantile of a sub-population
	,nSamplePop=getNSamples(mc)	##<< integer vector (nPop): number of samples in each population
	,nChainPop=getNChainsPop(mc)##<< number of chains within one population 
	, critSpecVarRatio=20		##<< if proprotion of spectral Density to Variation is higher than this value, signal problems and resort to subspaces 
	, dumpfileBasename="recover"##<< what to do on errors
	, pCheckSkipPart=0.5		##<< when checking each population for convergence problems, skip this proportion (kind of burnin) 
	, minNSampleCheck=16		##<< if a population has less samples within one chains, skip the check because it is too uncertain
){
	#nSamplePopChains <- nSamplePop*nChainPop  # number of samples across chains
	specVarRatioPop <- numeric(length(nSamplePop))	# initialized by 0
	#iPop=1
	for( iPop in seq_along(mc$pops)) if( nSamplePop[iPop] >= minNSampleCheck){
		# see .debug.divideTwDEMC below for commands to trace the problems
		pop <- if( pCheckSkipPart != 0 ){
			.subsetTwDEMCPop( mc$pops[[iPop]], iKeep= round((nSamplePop[iPop]-1)*pCheckSkipPart+1):nSamplePop[iPop] )			
		}else{
			mc$pops[[iPop]]
		}
		#mtrace(.checkProblemsSpectralPop)
		resFCheck <- .checkProblemsSpectralPop( pop, critSpecVarRatio= critSpecVarRatio)
#		if( resFCheck$hasProblems && (pSubs[iPop]/2 < minPSub) ){
#			if( !is.null(dumpfileBasename) )
#				if( dumpfileBasename == "recover"){
#					# see .debug.divideTwDEMC below for debuggin commands
#					cat("checkConvergenceProblems: encountered convergence problems of an already small population. calling recover() \n ")
#					recover()
#				}  else dump.frames(dumpfileBasename,TRUE)
#			stop(paste("checkConvergenceProblems: checking function encountered problems. Dump in ",dumpfileBasename,".rda",sep=""))
#		}
		specVarRatioPop[iPop] <- max(resFCheck$specVarRatioLumped)
	}
	### numeric vector (nPop): maximum ratios of spectral density to variance across parameters for each population
	specVarRatioPop
}

.debug.divideTwDEMC <- function(){
	# when encountering convergence problems in divideTwDMCSteps, recover is called
	# these functions then might be helpful
	mcp <- concatPops(subPops(mc, iPops=iPop)) 
	plot(as.mcmc.list(mcp), smooth=FALSE)
	(tmp <- .checkProblemsSpectralPop(pop,critSpecVarRatio= critSpecVarRatio))
	str3(mcp)
	getNGen(mcp) 
	mtrace(coda:::effectiveSize)
	effectiveSize( as.mcmc.list(mcp) )
	matplot( mcp$pAccept[,1,], type="l" )
	matplot( mcp$logDen[,1,], type="l" )
	matplot( mcp$temp, type="l" )
	ss <- stackChains(mcp)
	plot(ss[,"a"], ss[,"b"], col=rev(heat.colors(100))[ twRescale(ss[,1],c(20,100)) ] )
	#plot(mcp$Y[,"a",1], mcp$Y[,"b",1], col=rev(heat.colors(100))[ twRescale(mcp$Y[,"logDen1",1],c(20,100)) ] ) # -color does not scale to -Inf or very large values
	#plot(mcp$Y[,"a",1], mcp$Y[,"b",1] ) # -color does not scale to -Inf or very large values
	mcp$Y[,,1]
	#twCalcLogDenPar( den2dCor, mcp$Y[,1:2,1] ) # it is not -Inf but finite
	
	resFCheck
	
	tmp <- findSplit( ss[,-1], isBreakEarly=FALSE, rVarCrit = 2^2 )
	mtrace(getSubSpaces)
	tmp2 <- getSubSpaces( ss[,-1], isBreakEarly=FALSE, argsFSplit=list(rVarCrit = 2^2, rAlphaSlopeCrit=base:::pi/2),minPSub=ctrlSubspaces$minPSub, pSub=pSubs[iSub] )
	lapply(tmp2$spaces, "[[", "upperParBounds" )
	colnames(ss)[-1][tmp$iVars]
	
	iSubsCheckSubspaces <- which(  (pSubs > 2*ctrlSubspaces$minPSub) & (subPercChange > subPercChangeCrit) )
	#iSub <- iSubsCheck[length(iSubsCheck)]
	if( iSub %in% iSubsCheckSubspaces ){
		subSubSpaces <- getSubSpaces( ssp, isBreakEarly=FALSE, argsFSplit=argsFSplitPop[[ popSub]],minPSub=ctrlSubspaces$minPSub, pSub=pSubs[iSub] )
	}
	
	
	pSub > 2*ctrlSubspaces$minPSub
}

.splitPops <- function(
	### splitting subPopulations
	mcNewMinN				##<< twDEMCPops: last 2*m0 sample rows of last step
	, mcApp					##<< twDEMCPops: past and new samples
	, pSubs					##<< numeric vector (nPop): current percentile of the populations. Must sum to 1 within spaces
	, ctrlSubspaces			##<< list as in \code{\link{divideTwDEMCSACont}}
	, iPopsSplit				##<< integer vector (nPop): populations to split
	#, subPercChange			##<< numeric vector (nPop): relative change of subPopulations
	#, specVarRatioPop		##<< numeric vector (nPop): ratio of spectral denstity to varaince
	#, critSpecVarRatio=20	##<< if proprotion of spectral Density to Variation is higher than 0.8 times this value, split the population 
	#, subPercChangeCrit=1.6	##<< if all subPercChange of all sub-populations are below this value in two last batches, may assume convergence and skip further batches
	#, iPopsDoSplit=integer(0)	##<< populations given here will definitely be splitted (no matter of other criteria)
	#, minPSub = 0.1			##<< minimal quantile of a sub-population
	#, maxNSample=256		##<< if given a value, then looking for subspaces is done on a subsample of given length for efficiency (see \code{\link{getSubSpaces}})
	, nChainPop = getNChainsPop(mcApp)		##<< number of chains within population
	, m0 = calcM0twDEMC(getNParms(mcApp),nChainPop)	##<< minimum number of samples in step for extending runs
	, spaceInds = unique(getSpacesPop(mcApp))	##<< integer vector: set of space indices
	#, argsFSplit=list()	##<< further arguments to \code{\link{findSplit}}
){
	#print(".splitPop 1"); recover()
	pSubs1 <- pSubs		# 
	if( 0 != length(iPopsSplit) ){
		mcApp0 <- mcApp; mcNewMinN0 <- mcNewMinN	# save intial populations
		#mcApp <- mcApp0; mcNewMinN <- mcNewMinN0; pSubs1 <- pSubs
		nPop <- length(mcNewMinN$pops)
		nPopsSplit <- integer( nPop )	# number of subs that the population was splitted into
		#iPop = iCheckSplit[ length(iCheckSplit) ]
		newPops <- newPopsMcNew <- list()
		newPSubs <- integer(0)
		#iPop <- iPopsSplit[1]
		#iPop <- iPopsSplit[2]
		#iPop <- iPopsSplit[4]
		for( iPop in iPopsSplit){
			popMc <- mcNewMinN$pops[[iPop]]			# holds only the new samples
			aSample <- stackChains( popMc$parms )	# base splitting only on samples obtained in last batch
			#mtrace(getSubSpaces)
			subs <- getSubSpaces( aSample, isBreakEarly=FALSE, pSub=pSubs[iPop]
				#, minPSub=max(ctrlSubspaces$minPSub, pSubs[iPop]*(m0*nChainPop+(nChainPop-1))/min(nrow(aSample),ctrlSubspaces$maxNSample) )	# avoid populations with too few samples, regard loosing samples during splitting
				, minPSub=max(ctrlSubspaces$minPSub, pSubs[iPop]*(m0*nChainPop+(nChainPop-1))/min(nrow(aSample) ))  # maxNSample may be lower, however	
				, maxNSample=ctrlSubspaces$maxNSample
				, argsFSplit=ctrlSubspaces$argsFSplit	##<< further arguments to \code{\link{findSplit}}
				, splits=popMc$splits, lowerParBounds=popMc$lowerParBounds, upperParBounds=popMc$upperParBounds )
			# sapply( lapply( subs$spaces, "[[", "sample"), nrow )
			# ceiling(sapply( subs$spaces, "[[", "pSub")*maxNSample)
			#.getParBoundsPops(c(list(popMc),subs$spaces))
			if( length(subs$spaces) > 1){
				nPopsSplit[iPop] <- length(subs$spaces)	# mark splitted
				# split the overall samples
				newPSubs <- c( newPSubs, sapply(subs$spaces, "[[", "pSub") )
				pop <- mcApp$pops[[iPop]]
				#mtrace(divideTwDEMCPop)
				newPopsI <- divideTwDEMCPop(pop, subs$spaces)
				newPops <- c( newPops, newPopsI)
				# split the last 2*m0 rows sample
				popMcNew <- mcNewMinN$pops[[iPop]]
				newPopsMcNewI <- divideTwDEMCPop(popMcNew, subs$spaces)
				newPopsMcNew <- c( newPopsMcNew, newPopsMcNewI)
				#sapply( lapply( subs$spaces, "[[", "sample"), nrow )
				#sapply( lapply( newPopsI, "[[", "parms"), nrow )
				#sapply( lapply( newPopsMcNewI, "[[", "parms"), nrow )
				# check for consistency
				tmp <- lapply( newPopsI, .checkPop, nGen=12 )
				.nS <- nrow(pop$parms)
				.nSNew <- sum(sapply( newPopsI, function(jPop){ nrow(jPop$parms) }))
				if( (.nSNew > .nS) || (.nSNew < .nS-(nChainPop-1)*length(newPops)) )
					stop("splitPops: sample number of mcApp does not match after splitting.")
				# check for consistency
				tmp <- lapply( newPopsMcNewI, .checkPop, nGen=12 )
				.nS <- nrow(popMcNew$parms)
				.nSNew <- sum(sapply( newPopsMcNewI, function(jPop){ nrow(jPop$parms) }))
				if( (.nSNew > .nS) || (.nSNew < .nS-(nChainPop-1)*length(newPopsMcNewI)) )
					stop("splitPops: sample number of mcNewMinN does not match after splitting.")
				#.getParBoundsPops(c(newPopsI, list(pop)))
				#.getParBoundsPops(c(subs$spaces, list(pop))) #only internal splits
			}
		}	# for iSplit
		# actually replace the old by the new pops
		iPopsSplitted <- which(nPopsSplit != 0) 
		if( 0 != length(iPopsSplitted)){
			mcApp$pops <- c( mcApp0$pops[-iPopsSplitted], newPops )
			mcNewMinN$pops <- c( mcNewMinN$pops[-iPopsSplitted], newPopsMcNew )
			pSubs1 <- c( pSubs[-iPopsSplitted], newPSubs)
		}
		#
		# check consistency after splitting		
		.spacesPop <- getSpacesPop(mcApp)
		tmp <- lapply( mcApp$pops, .checkPop, nGen=12 )		
		if( !isTRUE( all.equal(.pSpaces <- as.numeric(tapply(pSubs1, .spacesPop, sum)), rep(1,length(unique(.spacesPop))))) )
			stop("splitPops (2c): pSubs within subspaces do not sum to 1.")
		#
		.nS1 <- sum(getNSamples(mcApp0))		# initial number of samples
		.nS2 <- sum(getNSamples(mcApp))
		#if( (.nS2 > .nS1) || (.nS2 < .nS1-(nChainPop-1)*length(iPopsSplit)) )
		if( (.nS2 > .nS1) || (.nS2 < .nS1-(nChainPop-1)*(getNPops(mcApp)-getNPops(mcApp0))) )
			stop("splitPops (2d): sample number does not match after splitting.")
		# seek neighboring populations
		for( iiSpace in 1:getNSpaces(mcApp) ){
			popsS <- mcApp$pops[ .spacesPop == spaceInds[iiSpace] ]
			.nPopS <- length(popsS)
			.iPopsS <- seq_along(popsS)
			if( .nPopS > 1) for( iPop in  .iPopsS){
					.popSource <- popsS[[iPop]]
					jPops <- .iPopsS[-iPop][ sapply( popsS[-iPop], function(pop){ 
								(length(pop$splits) >= length(.popSource$splits)) && all(pop$splits[ 1:length(.popSource$splits)] == .popSource$splits)
							})]
					if( length(jPops)==0 ) stop("splitPops: encountered no-neighboring case.")
					tmp.f <- function(){
						lapply(popsS, "[[", "splits")
						.spacePop <- getSpacesPop(mcApp0)
						lapply(subPops(mcApp0, iSpace=spaceInds[1])$pops, "[[", "splits")
						
					}
				}
		}
		#.getParBoundsPops(newPops)
	} # if length(iSplit)
	##value<< list with entries
	list(
		##describe<<
		mcApp = mcApp		##<< twDEMCPop with splitted populations
		,mcNewMinN = mcNewMinN	##<< twDEMCPop of last rows with splitted populations	
		,pSubs = pSubs1		##<< updated percentiles of the subPopulations. Sum to 1.
		)
		##end<<
}

twMergePops <- function(
	### Merge subspaces that hold only few samples or have only low percentiles to neighboring populations
	mcApp											##<< twDEMCPop: holding populations to merge
	#,mcNew=NULL	##<< second twDEMCPop: holding populations to merge
	, m0 = calcM0twDEMC(getNParms(mcApp),nChainPop)	##<< minimum number of samples 
	, minPSub = 0									##<< minimal quantile of a sub-population
	, pSubs	= getNSamples(mcApp)/getNSamplesSpace(mcApp)[getSpacesPop(mcApp)]	##<< numeric vector (nPops): percentiles of each population during last sampling (sum to 1) 	
	, nChainPop = getNChainsPop(mcApp)				##<< integer scalar: number of chains per population, maybe passed for efficiency reasons
	, spaceInds = unique(getSpacesPop(mcApp))		##<< integer vector: set of space indices, maybe passed for efficiency reasons
){
	mcApp0 <- mcApp; pSubs0 <- pSubs	# remeber state before merging
	#mcApp <- mcApp0;	pSubs <- pSubs0	# reset to initial for debugging
	fIPopsMerge <- function(mcApp,pSubs){	# make it a function to avoid replication and inconsitencies
		#which( pSubs < minPSub/2 | getNSamples(mcApp) < m0 )
		which( getNSamples(mcApp) < m0 )
	}
	iPopsMerge <- iPopsMerge0 <- fIPopsMerge(mcApp,pSubs)
	iExitMerge <- getNPops(mcApp0)	# prevent infinite runs
	.nSplits <- 0
	#if( 0 != length(iPopsMerge) ) recover()
	nSpaces <- getNSpaces(mcApp)
	iPopsBefore <- seq_along(mcApp$pops)		# at the start the indices correspond
	while( 0 != length(iPopsMerge) && iExitMerge != 0){
		iPop <- iPopsMerge[ order(pSubs[iPopsMerge])[1] ] # the one with lowest p
		#print(iPop)
		#str3(mcApp$pops[[iPop]])
		#mcApp$pops[[iPop]]$splits
		#iPopsSpace <- which(getSpacesPop(mcApp) == mcApp$pops[[iPop]]$spaceInd)
		#lapply(mcApp$pops[iPopsSpace], "[[", "splits")
		#lapply(mcApp$pops[iPopsSpace], "[[", "splits")
		#mtrace(.mergePopTwDEMC)
		mcAppPrev <- mcApp; pSubsPrev <- pSubs; iPopsBeforePrev <- iPopsBefore		#remember state before current merge
		#mcApp <- mcAppPrev; pSubs1 <- pSubsPrev
		#mtrace(.mergePopTwDEMC)
		resMerge <- .mergePopTwDEMC( mcApp$pops, iPop, pSubs, mergeMethod="random" ) # use random to avoid high spectral density calculation afterwards
		#resMergeNew <- .mergePopTwDEMC( mcNew$pops, iPop, pSubs, mergeMethod="random" ) # use random to avoid high spectral density calculation afterwards
		# consitency check: seek neighboring populations
		for( iSpace in 1:nSpaces ){
			popsS <- resMerge$pops[ sapply(resMerge$pops, "[[","spaceInd") == spaceInds[iSpace] ]
			.nPopS <- length(popsS)
			.iPopsS <- seq_along(popsS)
			if( .nPopS > 1) for( iPop in  .iPopsS){
					.popSource <- popsS[[iPop]]
					jPops <- .iPopsS[-iPop][ sapply( popsS[-iPop], function(pop){ 
								(length(pop$splits) >= length(.popSource$splits)) && all(pop$splits[ 1:length(.popSource$splits)] == .popSource$splits)
							})]
					if( length(jPops)==0 ) stop("mergePops_resMerge: encountered no-neighboring case.")
				}
		}
		#mcApp$pops[[iPop]]$spaceInd
		#rms <- sapply(resMerge$pops,"[[","spaceInd");	iPopsSpaceNew <- which(rms == 1); (tmp1<-lapply(resMerge$pops[iPopsSpaceNew], "[[", "splits"))
		#rms <- sapply(resMerge$pops,"[[","spaceInd");	iPopsSpaceNew <- which(rms == 2); (tmp2<-lapply(resMerge$pops[iPopsSpaceNew], "[[", "splits"))
		#print(sum(sapply(lapply(mcApp$pops,"[[","parms"),nrow)));	print(sum(sapply(lapply(resMerge$pops,"[[","parms"),nrow)))
		#
		#if( sum( sapply(lapply(resMerge$pops,"[[","parms"),nrow) ) > .nS2 ) { print("divideTwDEMCSteps: encountered larger sample size"); recover() }
		#
		#-- replace the pops
		mcApp$pops <- resMerge$pops
		#mcNew$pops <- resMergeNew$pops
		pSubs <- resMerge$pPops
		iPopsBefore <- iPopsBeforePrev[ resMerge$iPopsBefore ]	#update the mapping of new to old indices
		.nSplits <- .nSplits + length(resMerge$pPopsSource)
		iExitMerge <- iExitMerge -1 
		#
		#-- reevaluate the populations to merge
		iPopsMerge <- fIPopsMerge(mcApp,pSubs)
		#
		#-- check consistency
		# there may be more populations without samples, check only the ones that have samples
		tmp <- lapply( mcApp$pops[getNSamples(mcApp)!=0], .checkPop, nGen=12 )		
		.spacesPop <- getSpacesPop(mcApp)
		if( !isTRUE( all.equal(.pSpaces <- as.numeric(tapply(pSubs, .spacesPop, sum)), rep(1,length(unique(.spacesPop))))) )
			stop("mergePops (3): pSubs within subspaces do not sum to 1.")
		# seek neighboring populations
		for( iSpace in 1:getNSpaces(mcApp) ){
			popsS <- mcApp$pops[ .spacesPop == spaceInds[iSpace] ]
			.nPopS <- length(popsS)
			.iPopsS <- seq_along(popsS)
			if( .nPopS > 1) for( iPop in  .iPopsS){
					.popSource <- popsS[[iPop]]
					jPops <- .iPopsS[-iPop][ sapply( popsS[-iPop], function(pop){ 
								(length(pop$splits) >= length(.popSource$splits)) && all(pop$splits[ 1:length(.popSource$splits)] == .popSource$splits)
							})]
					if( length(jPops)==0 ) stop("mergePops_merge: encountered no-neighboring case.")
					tmp.f <- function(){
						lapply(popsS, "[[", "splits")
						.spacePop <- getSpacesPop(mcApp0)
						lapply(subPops(mcApp0, iSpace=spaceInds[1])$pops, "[[", "splits")
					}
				}
		}
		#
		#print(iPopsMerge)
	}# while
	if( iExitMerge == 0 ) stop("mergePops: while loop of populations to merge did not exit.")
	# check for consistency
	.nS1 <- sum(getNSamples(mcApp0))
	.nS3 <- sum(sapply( mcApp$pops, function(jPop){ nrow(jPop$parms) }))
	if( (.nS3 > .nS1) || (.nS3 < .nS1-(nChainPop-1)*.nSplits) )
		stop("mergePops (4): sample number does not match after merging.")
	if( any(getNSamples(mcApp) < m0) ){ print("twMergePops: Too few samples, recover"); recover() }
	##value<< list with entries
	list(
		##describe<<
		mcApp = mcApp		##<< twDEMCPop with merged populations
		#,mcNew = mcNew		##<< twDEMCPop with merged populations
		,pSubs = pSubs		##<< updated percentiles of the subPopulations. Sum to 1.
		,iPopsBefore = iPopsBefore	##<< integer vector (nPop): indices of the largest original population contributing to the new population 
	)
	##end<<
}








	





