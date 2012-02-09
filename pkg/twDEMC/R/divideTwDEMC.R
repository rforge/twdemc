# see also subspace.R

divideTwDEMCStep <- function(
	### run twDEMCBlock on subspaces
	aTwDEMC					##<< former run of twDEMCBlockInt
	,qPop=numeric(0)		##<< numeric vector (nPop): probability of samples in subspace within entire space 
	,nGen=512				##<< the number of generations for twDEMCBlock
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMCBlock}} containing entry thin
	, nSampleMin=16			##<< minimum number of samples in each population so that calculation of average density is stable
	, m0 = calcM0twDEMC(getNParms(aTwDEMC),getNChainsPop(aTwDEMC))		# minimum number of samples in step for extending runs
	, ...					##<< further arguments to \code{\link{twDEMCBlock}}
){
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
	nGen0PopsThin <- pmax( m0, nSampleMin, ceiling(nGen0Pops/thin))*thin	# must have at least m0 samples in order to be able to extend the sample	 
	#
	#----- initial run based on current quantiles
	#trace("twDEMCBlockInt")
	resTwDEMC <- resTwDEMC0 <- twDEMCBlock(aTwDEMC, nGen=nGen0PopsThin, extendRun=FALSE, controlTwDEMC=controlTwDEMC, m0=m0, ...)
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
	popLogMeanDensSubs <- lapply( 1:nSpace, function(iiSpace){
			iPops <- iPopsSpace[[iiSpace]]
			#iSub <- 2
			#iSub <- nSub
			lw <- sapply(iPops,function(iPop){  # average the unnormalized densities
					ss <- stackChains(resTwDEMC$pops[[iPop]]$logDen)
					ssLogDen <- rowSums(ss)
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
		p2u <- exp(popLogMeanDensSubs[[iiSpace]])*qPopSpace[[iiSpace]]		# two weights: quantile and density importance
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
			# stay at current temperature (we are alreay at the prescribed end temperature)
			for( iPop in seq_along(resTwDEMC0$pops) ){
				resTwDEMC0$pops[[iPop]]$TEnd <- 
					resTwDEMC0$pops[[iPop]]$temp[ nSamples0[iPop] ]
			}
			#mtrace(twDEMCBlock.twDEMCPops)
			#mtrace(twDEMCBlockInt)
			resTwDEMC <- twDEMCBlock( resTwDEMC0, nGen=nSamplesAdd*thin, controlTwDEMC=controlTwDEMC,  m0=m0, ...  )
		}
	#all( getNSamples(resTwDEMC) >= nSamplesSubsReq)
	#getNSamples(resTwDEMC1)
	
	#----- do a subsampling
	# because we required minNSamples to take we possibly can add more than required
	resTwDEMC <- resTwDEMC2 <- squeeze.twDEMCPops(resTwDEMC1, length.out=nSamplesSubsReq)
	#getNSamples(resTwDEMC2)

	#----- append to former run
	#resTwDEMC <- .concatTwDEMCRuns(aTwDEMC,resTwDEMC2,doRecordProposals=doRecordProposals)
	#getNSamples(resTwDEMC)

	#------ check for consistency
	#lapply( )
	
	##value<< 
	## For each population, a list with entries
	res <- list(	##describe<<
		resTwDEMC = resTwDEMC	##<< the runs with new samples of class twDEMCPops
		,pSubs = pPops			##<< the quantiles of samples within given subspace
		,subPercChange = subPercChange	##<< ##<< numeric vector: ratio of estimated proportion in limiting distribution to proportion of initial proportion
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
		specVarRatioLumped <- spec/varSample
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
				spec/varSample
			})
		list(
			hasProblems = TRUE
			,specVarRatioLumped = specVarRatioLumped
			,specVarRatio = specVarRatio 
		)
	}
}

divideTwDEMCSACont <- function(
	### run twDEMCBlock on subspaces with decreasing Temperature
	mc					##<< former run of twDEMCBlockInt
	, nObs				##<< integer vector (nResComp) specifying the number of observations for each result component
	, nGen=512			##<< the number of generations for twDEMCBlock
	#, nBatch=4			##<< number of batches with recalculated Temperature
	, ...				##<< further arguments to \code{\link{twDEMCBlock}}
	, TFix=numeric(0)	##<< numeric vector (nResComp) specifying a finite value for components with fixed Temperatue, and a non-finite for others
	, TMax=numeric(0)	##<< numeric vector (nResComp) specifying a maximum temperature for result components.
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMCBlock}} containing entry thin
	, doRecordProposals=FALSE	##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y
	, m0 = calcM0twDEMC(getNParms(mc),getNChainsPop(mc))	##<< minimum number of samples in step for extending runs
	, nGenBatch=512											##<< number of runs for one twDEMCStep
	#, argsFSplitPop=vector("list",dim(aSample)[3])			##<< for each population: list of arguments  passed \code{\link{getSubSpaces}} and further to \code{\link{findSplit}}, e.g. for passing order of variables to check in \code{iVars} and \code{jVarsVar}
	#
	, dumpfileBasename="recover"##<< what to do on errors
	, minPSub = 0.1				##<< minimal quantile of a sub-population
	, subPercChangeCrit=1.6		##<< if all subPercChange of all sub-populations are below this value in two last batches, may assume convergence and skip further batches
	, maxNSample=256			##<< if given a value, then looking for subspaces is done on a subsample of given length for efficiency (see \code{\link{getSubSpaces}})
	, pThinPast=0.5				##<< in each step thin the past to given fraction before appending it
	, maxRelTChange=0.025 		##<< if Temperature of the components changes less than specified value, the algorithm can finish
	, maxLogDenDrift=0.3		##<< if difference between mean logDensity of first and fourth quartile of the sample is less than this value, we do not need further batches because of drift in logDensity
	, TDecProp=0.9				##<< proportion of Temperature decrease: below one to diminish risk of decreasing Temperature too fast (below what is supported by other data streams)
	, gelmanCrit=1.4			##<< do not decrease Temperature, if variance between chains is too high, i.e. Gelman Diag is above this value
	, critSpecVarRatio=20		##<< if proprotion of spectral Density to Variation is higher than this value, signal problems and resort to subspaces
	, pCheckSkipPart=0.5		##<< when checking each population for convergence problems, skip this proportion (kind of burnin) 
	, minNSampleCheck=32		##<< if a population has less samples across chains, skip the check because it is too uncertaint
	, restartFilename=NULL		##<< filename to write intermediate results to
	, iPopsDoSplit=integer(0)	##<< populations given here will definitely be splitted (no matter of other criteria)
){
	args <- c( list( nObs=nObs, TFix=TFix, maxRelTChange=maxRelTChange, maxLogDenDrift=maxLogDenDrift,  TDecProp=TDecProp, restartFilename=restartFilename)
		, list(controlTwDEMC=controlTwDEMC, doRecordProposals=doRecordProposals, m0=m0, nGenBatch=nGenBatch) 
		, list(...) 
	)
	nResComp <- ncol(mc$pops[[1]]$resLogDen)
	if( 0 == length( TFix) ) TFix <- structure( rep( NA_real_, nResComp), names=colnames(mc$pops[[1]]$resLogDen) )
	if( nResComp != length( TFix) ) stop("divideTwDEMCSACont: TFix must be of the same length as number of result Components.")
	iFixTemp <- which( is.finite(TFix) )
	if( 0 == length( TMax) ) TMax <- structure( rep( NA_real_, nResComp), names=colnames(mc$pops[[1]]$resLogDen) )
	if( nResComp != length( TMax) ) stop("divideTwDEMCSACont: TMax must be of the same length as number of result Components.")
	iMaxTemp <- which( is.finite(TMax) )
	#
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
	thin <- controlTwDEMC$thin
	if( nGenBatch%/%thin < 64) 
		stop("not enough sample per batch. Increase nGenBatch to at least 64*thin or decreaste ctrolTwDEMC$thin")
	nGenBatch = (nGenBatch%/%thin)*thin		# make nGen multiple of thin
	boPopsDidConverge <- FALSE
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
	spacePop <- getSpacesPop(mc)	# integer vector (nPop): specifying the space replicated that the population belongs to
	spaceInds <- unique(spacePop)
	nSamplesPop <- getNSamples(mc)
	nSamplesSpace <- structure( sapply(spaceInds, function(spaceInd){ sum(nSamplesPop[spacePop==spaceInd])}), names=spaceInds)
	#
	# the following updated variables are required in each cycle
	mcNew <- mc		# new samples within batch
	mcApp <- mc		# concatenation of thinned past and new samples
	pSubs <- nSamplesPop/nSamplesSpace[ as.character(spacePop) ]
	subPercChange=numeric(length(pSubs))	# initialized to zero
	specVarRatioPop <- .checkConvergenceProblems(mcApp, critSpecVarRatio=Inf, pCheckSkipPart=pCheckSkipPart)
	#
	while( (iGen < nGen) && !boPopsDidConverge){
		iBatch <- iGen/nGenBatch+1
		mcApp0 <- mcApp; pSubs0 <- pSubs	# remember the initial populations
		#mcApp<-mcApp0; pSubs<-pSubs0
		.saveRestartFile( restartFilename, mcApp=mcApp, args=args )
		#
		#-- can decrease Temperature be decreased?
		TCurr <- getCurrentTemp(mcApp) 
		mcSpace <- stackPopsInSpace(mcApp, mergeMethod="random")	
		mcSpaceEnd <- thin(mcSpace, start=getNGen(mcSpace) %/% 2)
		mcl <- as.mcmc.list(stackChainsPop( mcSpaceEnd ))
		#print("divideTwDEMCSACont: loop start"); recover()
		#plot( mcl, smooth=FALSE )
		#plot( res$pops[[1]]$parms[,"b",1] ~ res$pops[[1]]$parms[,"a",1], col=rainbow(255)[round(twRescale(-res$pops[[1]]$resLogDen[,"logDen1",1],c(10,200)))] )
		gelmanDiagRes <- try( gelman.diag(mcl)$mpsrf )	# cholesky decomposition may throw errors
		print(paste("gelmanDiag=",signif(gelmanDiagRes,2)," T=",paste(signif(TCurr,2),collapse=","), sep="") )
		TEnd <- .calcTEnd(gelmanDiagRes=gelmanDiagRes, resEnd=mcSpaceEnd, TCurr=TCurr, nObsDen=nObsDen, TFix=TFix, iFixTemp=iFixTemp, iNonFixTempDens=iNonFixTempDens, TMax=TMax, iMaxTemp=iMaxTemp, gelmanCrit=gelmanCrit, TDecProp=TDecProp)
		#
		#-- split pops that changed poportion a lot or pops with high spectralDensity ratio into smaller pops
		#trace(.splitPops, recover )
		resSplitPops <- .splitPops(mcNew=mcNew, mcApp=mcApp, pSubs, subPercChange=subPercChange, specVarRatioPop=specVarRatioPop
			,minPSub=minPSub, subPercChangeCrit=subPercChangeCrit, maxNSample=maxNSample, m0=m0, nChainPop=nChainPop, spaceInds=spaceInds
			,iPopsDoSplit=iPopsDoSplit
		)
		mcApp <- mcApp1 <- resSplitPops$mcApp
		pSubs <- pSubs1 <- resSplitPops$pSubs
		iPopsDoSplit <- integer(0)
		#
		#-- sample the next batch
		nGenStep <- min(nGenBatch, nGen-iGen)
		#mtrace(divideTwDEMCStep)
		resStep <- divideTwDEMCStep(mcApp, nGen=nGenStep, doRecordProposals=doRecordProposals, m0=m0, TEnd=TEnd, ... )
		mcNew <- resStep$resTwDEMC			# mcNew: only the new samples 
		pSubs <- pSubs2 <- pSubsNew <- resStep$pSubs	
		subPercChange <- resStep$subPercChange
		.spacesPop <- getSpacesPop(mcNew)
		if( !isTRUE( all.equal(.pSpaces <- as.numeric(tapply(pSubs, .spacesPop, sum)), rep(1,length(unique((.spacesPop)))) )) )
			stop("divideTwDEMCSteps (1): pSubs within subspaces do not sum to 1.")
		#
		#-- check convergence problems of the new samples
		#specVarRatioPop <- .checkConvergenceProblems(mcNew, critSpecVarRatio=Inf, pCheckSkipPart=pCheckSkipPart)
		specVarRatioPop <- .checkConvergenceProblems(mcNew, nChainPop=nChainPop, critSpecVarRatio=critSpecVarRatio, dumpfileBasename=dumpfileBasename, pCheckSkipPart=pCheckSkipPart, minNSampleCheck=minNSampleCheck)
		#
		#-- append to thinned past
		#.tmp <- lapply( mcNew$pops[nSamplePop != 0], .checkPop, nGen=12 )		
		#.spacesPop <- getSpacesPop(mcNew)
		#tapply(nSamplePop, .spacesPop, function(nsI){ nsI/sum(nsI)})
		#
		#.spacesPopOrig <- getSpacesPop(mcApp)
		#.pSpacesOrig <- as.numeric(tapply(pSubs, .spacesPopOrig, sum)
		#boPopsDistConverge <- all(resStep$subPercChange <= subPercChangeCrit)
		mcAppPast <- squeeze(mcApp, length.out=ceiling(getNSamples(mcApp)*pThinPast) ) # mcApp0: thinned former mcApp
		mcApp <- mcApp2 <- .concatTwDEMCRuns(mcAppPast,mcNew,doRecordProposals=doRecordProposals) # mcApp and mcApp1: thinned past + new samples 
		#.nS1 <- sum(getNSamples(mcApp2))
		#if( nGenStep/thin >= 64){  # only check and split populations, if enough samples in batch
		#iPop=length(mcNew$pops)
		#
		#-- merge subspaces that contain less samples than minPSub/2
		# again base decision on sample of last batch (pSubs1)
		resMergePops <- .mergePops(mcApp=mcApp, pSubs=pSubs, minPSub=minPSub, nChainPop=nChainPop, m0=m0, spaceInds=spaceInds)
		mcApp <- mcApp3 <- resMergePops$mcApp
		pSubs <- pSubs3 <- resMergePops$pSubs 
		if( any(getNSamples(mcApp) < m0) ){ print("divideTwDEMCSACont: Too few samples, recover"); recover() }
		#} # enough new samples for checking/splitting
		#
		iGen <- iGen + nGenStep
		#cat(paste(iGen," out of ",nGen," generations completed. T=",paste({T<-res$temp[nrow(res$temp),];round(T,digits=ifelse(T>20,0,1))},collapse="  "),"     ",date(),"\n",sep=""))
		cat(paste(iGen," out of ",nGen," generations completed. ,max(subPercChange)=",max(resStep$subPercChange),"     ",date(),"\n",sep=""))
	}	# while iGen < nGenBatch
	resStep$resTwDEMC <- mcApp
	resStep
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
	#trace(divideTwDEMCSACont, recover)
	resDiv <- res2Div <- divideTwDEMCSACont( res, nObs=nObs, nGen=1024 
		, debugSequential=TRUE
		, controlTwDEMC=list(DRgamma=0.1)		# DR step of 1/10 of the proposed length
		, iPopsDoSplit=1:2		# definitely split the populations
	)
	resDiv <- res3Div <- divideTwDEMCSACont( resDiv$resTwDEMC, nObs=nObs, nGen=2048 )
	resDiv$pSubs
	resDiv$subPercChange
	getNSamples(resDiv$resTwDEMC)
	res2 <- stackChainsPop(tmp <- stackPopsInSpace(resDiv$resTwDEMC, mergeMethod="random", nInSlice = 1))
	plot( tmp$pops[[1]]$temp )
	#plot( res2$pops[[1]]$temp )
	getCurrentTemp(res2)
	plot( as.mcmc.list(res2), smooth=FALSE )
	mc <- concatPops(res2)
	matplot( mc$pAccept[,1,], type="l" )
	matplot( mc$resLogDen[,1,], type="l" )
	iSpace=2; plot( mc$parms[,"a",iSpace], mc$parms[,"b",iSpace], xlim=c(-0.5,2), ylim=c(-20,20), col=(heat.colors(100))[twRescale(log(-mc$resLogDen[,1,iSpace]),c(10,100))])
	gelman.diag(as.mcmc.list(res2))
	.checkProblemsSpectralPop(res2$pops[[1]])
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
	TEnd <- if(  !inherits(gelmanDiagRes,"try-error") && gelmanDiagRes <= gelmanCrit ){
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
			TEnd <- TCurr - TDecProp*(TCurr-TEnd)
		}else{
			TEnd <-TCurr
		}
	### numeric vector (nResComp): target temperature for all result components 
	TEnd
}

.checkConvergenceProblems <- function(
	### checking all population for convergence problems
	mc			##<< result of twDEMC, usually only the last samples without the history
	,nSamplePop=getNSamples(mc)	##<< integer vector (nPop): number of samples in each population
	,nChainPop=getNChainsPop(mc)##<< number of chains within one population 
	, critSpecVarRatio=20		##<< if proprotion of spectral Density to Variation is higher than this value, signal problems and resort to subspaces 
	, dumpfileBasename="recover"##<< what to do on errors
	, pCheckSkipPart=0.5		##<< when checking each population for convergence problems, skip this proportion (kind of burnin) 
	, minNSampleCheck=32		##<< if a population has less samples across chains, skip the check because it is too uncertaint
){
	nSamplePopChains <- nSamplePop*nChainPop  # number of samples across chains
	specVarRatioPop <- numeric(length(nSamplePop))	# initialized by 0
	for( iPop in seq_along(mc$pops)) if( nSamplePopChains[iPop] >= minNSampleCheck){
		# see .debug.divideTwDEMC below for commands to trace the problems
		pop <- if( pCheckSkipPart != 0 ){
			.subsetTwDEMCPop( mc$pops[[iPop]], iKeep= round((nSamplePop[iPop]-1)*pCheckSkipPart+1):nSamplePop[iPop] )			
		}else{
			mc$pops[[iPop]]
		}
		#mtrace(.checkProblemsSpectralPop)
		resFCheck <- .checkProblemsSpectralPop( pop, critSpecVarRatio= critSpecVarRatio)
		if( resFCheck$hasProblems ){
			if( !is.null(dumpfileBasename) )
				if( dumpfileBasename == "recover"){
					# see .debug.divideTwDEMC below for debuggin commands
					cat("checkConvergenceProblems: encountered convergence problems in a subSpace. calling recover() \n ")
					recover()
				}  else dump.frames(dumpfileBasename,TRUE)
			stop(paste("checkConvergenceProblems: checking function encountered problems. Dump in ",dumpfileBasename,sep=""))
		}
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
	(tmp <- .checkProblemsSpectralPop(pop))
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
	tmp2 <- getSubSpaces( ss[,-1], isBreakEarly=FALSE, argsFSplit=list(rVarCrit = 2^2, rAlphaSlopeCrit=base:::pi/2),minPSub=minPSub, pSub=pSubs[iSub] )
	lapply(tmp2$spaces, "[[", "upperParBounds" )
	colnames(ss)[-1][tmp$iVars]
	
	iSubsCheckSubspaces <- which(  (pSubs > 2*minPSub) & (subPercChange > subPercChangeCrit) )
	#iSub <- iSubsCheck[length(iSubsCheck)]
	if( iSub %in% iSubsCheckSubspaces ){
		subSubSpaces <- getSubSpaces( ssp, isBreakEarly=FALSE, argsFSplit=argsFSplitPop[[ popSub]],minPSub=minPSub, pSub=pSubs[iSub] )
	}
	
	
	pSub > 2*minPSub
}

.splitPops <- function(
	### splitting subPopulations
	mcNew					##<< twDEMCPops: new samples of last step
	, mcApp					##<< twDEMCPops: past and new samples
	, pSubs					##<< numeric vector (nPop): current percentile of the populations. Must sum to 1 within spaces
	, subPercChange			##<< numeric vector (nPop): relative change of subPopulations
	, specVarRatioPop		##<< numeric vector (nPop): ratio of spectral denstity to varaince
	, critSpecVarRatio=20	##<< if proprotion of spectral Density to Variation is higher than 0.8 times this value, split the population 
	, minPSub = 0.1			##<< minimal quantile of a sub-population
	, subPercChangeCrit=1.6	##<< if all subPercChange of all sub-populations are below this value in two last batches, may assume convergence and skip further batches
	, maxNSample=256		##<< if given a value, then looking for subspaces is done on a subsample of given length for efficiency (see \code{\link{getSubSpaces}})
	, nChainPop = getNChainsPop(mcNew)		##<< number of chains within population
	, m0 = calcM0twDEMC(getNParms(mcNew),nChainPop)	##<< minimum number of samples in step for extending runs
	, spaceInds = unique(getSpacesPop(mcApp))	##<< integer vector: set of space indices
	, iPopsDoSplit=integer(0)	##<< populations given here will definitely be splitted (no matter of other criteria)
){
	pSubs1 <- pSubs		# 
	iPopsSplit <- union( iPopsDoSplit, which( (pSubs/2 > minPSub) & ((subPercChange >= subPercChangeCrit) | (specVarRatioPop > 0.8*critSpecVarRatio)) ))
	if( 0 != length(iPopsSplit) ){
		mcApp0 <- mcApp	# save intial populations
		#iPop = iCheckSplit[ length(iCheckSplit) ]
		newPops <- list()
		newPSubs <- integer(0)
		#iPop <- iPopsSplit[1]
		for( iPop in iPopsSplit){
			popMc <- mcNew$pops[[iPop]]			# mc holds only the new samples
			aSample <- stackChains( popMc$parms )	# base splitting only on samples obtained in last batch
			#mtrace(getSubSpaces)
			subs <- getSubSpaces( aSample, isBreakEarly=FALSE, pSub=pSubs[iPop]
				, minPSub=max(minPSub, (m0*nChainPop+(nChainPop-1))/maxNSample )	# avoid populations with too few samples, regard loosing samples during splitting
				, maxNSample=maxNSample
				, splits=popMc$splits, lowerParBounds=popMc$lowerParBounds, upperParBounds=popMc$upperParBounds )
			# sapply( lapply( subs$spaces, "[[", "sample"), nrow )
			# ceiling(sapply( subs$spaces, "[[", "pSub")*maxNSample)
			#.getParBoundsPops(c(list(popMc),subs$spaces))
			if( length(subs$spaces) > 1){
				pop <- mcApp$pops[[iPop]]
				#mtrace(divideTwDEMCPop)
				newPopsI <- divideTwDEMCPop(pop, subs$spaces)
				#sapply( lapply( newPopsI, "[[", "parms"), nrow )
				newPops <- c( newPops, newPopsI)
				newPSubs <- c( newPSubs, sapply(subs$spaces, "[[", "pSub") )
				# check for consistency
				lapply( newPopsI, .checkPop, nGen=12 )
				.nS <- nrow(pop$parms)
				.nSNew <- sum(sapply( newPopsI, function(jPop){ nrow(jPop$parms) }))
				if( (.nSNew > .nS) || (.nSNew < .nS-(nChainPop-1)*length(newPops)) )
					stop("splitPops: sample number does not match after splitting.")
				#.getParBoundsPops(c(newPopsI, list(pop)))
				#.getParBoundsPops(c(subs$spaces, list(pop))) #only internal splits
			}else{
				iPopsSplit[match(iPop,iPopsSplit)] <- NA		# remove from pops to split (do not remove in mcApp)
			}
		}
		# actually replace the old by the new pops
		iPopsSplit <- iPopsSplit[ is.finite(iPopsSplit) ] 
		if( 0 != length(iPopsSplit)){
			mcApp$pops <- c( mcApp0$pops[-iPopsSplit], newPops )
			pSubs1 <- c( pSubs[-iPopsSplit], newPSubs)
		}
		#
		# check consistency after splitting		
		.spacesPop <- getSpacesPop(mcApp)
		lapply( mcApp$pops, .checkPop, nGen=12 )		
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
		,pSubs = pSubs1		##<< updated percentiles of the subPopulations. Sum to 1.
		)
		##end<<
}

.mergePops <- function(
	### merge populations with very low share
	mcApp		##<< twDEMCPop: holding populations to merge
	,pSubs		##<< numeric vector (nPops): percentiles of each population during last sampling (sum to 1) 	
	, minPSub = 0.1			##<< minimal quantile of a sub-population
	, nChainPop = getNChainsPop(mcApp)
	, m0 = calcM0twDEMC(getNParms(mcApp),nChainPop)	##<< minimum number of samples in step for extending runs
	, spaceInds = unique(getSpacesPop(mcApp))	##<< integer vector: set of space indices
){
	mcApp0 <- mcApp; pSubs0 <- pSubs	# remeber state before merging
	#mcApp <- mcApp0;	pSubs1 <- pSubs0	# reset to initial for debugging
	iPopsMerge <- iPopsMerge0 <- which( pSubs < minPSub/2 | getNSamples(mcApp) < m0 )
	iExitMerge <- getNPops(mcApp0)	# prevent infinite runs
	.nMerges <- 0
	#if( 0 != length(iPopsMerge) ) recover()
	nSpaces <- getNSpaces(mcApp)
	while( 0 != length(iPopsMerge) && iExitMerge != 0){
		#print(iPop)
		iPop <- iPopsMerge[ order(pSubs[iPopsMerge])[1] ] # the one with lowest p
		#str3(mcApp$pops[[iPop]])
		#mcApp$pops[[iPop]]$splits
		#iPopsSpace <- which(getSpacesPop(mcApp) == mcApp$pops[[iPop]]$spaceInd)
		#lapply(mcApp$pops[iPopsSpace], "[[", "splits")
		#lapply(mcApp$pops[iPopsSpace], "[[", "splits")
		#mtrace(.mergePopTwDEMC)
		mcAppPrev <- mcApp; pSubsPrev <- pSubs		#remember state before current merge
		#mcApp <- mcAppPrev; pSubs1 <- pSubsPrev
		#mtrace(.mergePopTwDEMC)
		resMerge <- .mergePopTwDEMC( mcApp$pops, iPop, pSubs, mergeMethod="random" ) # use random to avoid high spectral density calculation afterwards
		# seek neighboring populations
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
		pSubs <- resMerge$pPops
		.nMerges <- .nMerges + resMerge$nSubs
		iExitMerge <- iExitMerge -1 
		#
		#-- reevaluate the populations to merge
		iPopsMerge <- which( pSubs < minPSub/2 | getNSamples(mcApp) < m0)
		#
		#-- check consistency
		tmp <- lapply( mcApp$pops, .checkPop, nGen=12 )		
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
	}
	if( iExitMerge == 0 ) stop("mergePops: while loop of populations to merge did not exit.")
	# check for consistency
	.nS1 <- sum(getNSamples(mcApp0))
	.nS3 <- sum(sapply( mcApp$pops, function(jPop){ nrow(jPop$parms) }))
	if( (.nS3 > .nS1) || (.nS3 < .nS1-(nChainPop-1)*.nMerges) )
		stop("mergePops (4): sample number does not match after merging.")
	if( any(getNSamples(mcApp) < m0) ){ print(".mergePops: Too few samples, recover"); recover() }
	##value<< list with entries
	list(
		##describe<<
		mcApp = mcApp		##<< twDEMCPop with merged populations
		,pSubs = pSubs		##<< updated percentiles of the subPopulations. Sum to 1.
		,popMappint = TODO	##<< integer vector (nPop_old): giving the index of the first new population
	)
	##end<<
}








	





