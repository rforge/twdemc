divideTwDEMCSteps <- function(
	### run twDEMCBlock on subspaces
	aTwDEMC					##<< former run of twDEMCBlockInt
	,nGen=512				##<< the number of generations for twDEMCBlock
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMCBlockInt}} containing entry thin
	, doRecordProposals=FALSE		##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y
	, m0 = calcM0twDEMC(getNParms(aTwDEMC),getNChainsPop(aTwDEMC))		# minimum number of samples in step for extending runs
	, ...					##<< further arguments to \code{\link{twDEMCBlockInt}}
	, nGenBatch=256
	#, argsFSplitPop=vector("list",dim(aSample)[3])	##<< for each population: list of arguments  passed \code{\link{getSubSpaces}} and further to \code{\link{findSplit}}, e.g. for passing order of variables to check in \code{iVars} and \code{jVarsVar}
	, dumpfileBasename="recover"
	, critSpecVarRatio=25
	, minPSub = 0.1
	, subPercChangeCrit=1.6	##<< if all subPercChange of all sub-populations are below this value in two last batches, may assume convergence and skip further batches
	, maxNSample=256					##<< if given a value, then looking for subspaces is done on a subsample of given length for efficiency (see \code{\link{getSubSpaces}})
	, pThinPast=0.5		##<< in each step thin the past to given fraction before appending it
){
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
	thin <- controlTwDEMC$thin
	if( nGenBatch%/%thin < 64) 
		stop("not enough sample per batch. Increase nGenBatch to at least 64*thin or decreaste ctrolTwDEMC$thin")
	nGenBatch = (nGenBatch%/%thin)*thin		# make nGen multiple of thin
	boPopsDidConverge <- FALSE
	mcApp <- aTwDEMC		#new samples are appended each batch
	iGen = 0
	nChainPop <- getNChainsPop(aTwDEMC)
	while( (iGen < nGen) && !boPopsDidConverge){
		#------- sample the next batch
		nGenStep <- min(nGenBatch, nGen-iGen)
		#mtrace(divideTwDEMCStep)
		resStep <- divideTwDEMCStep(mcApp, nGen=nGenStep, doRecordProposals=doRecordProposals, m0=m0, ... )
		mc <- resStep$resTwDEMC			# mc: only the new samples 
		pSubs <- pSubs1 <- resStep$pSubs	#pSubs1 will reflect splitted populations
		#
		nSamplePop <- getNSamples(mc)
		.tmp <- lapply( mc$pops[nSamplePop != 0], .checkPop, nGen=12 )		
		.spacesPop <- getSpacesPop(mc)
		#tapply(nSamplePop, .spacesPop, function(nsI){ nsI/sum(nsI)})
		if( !isTRUE( all.equal(.pSpaces <- as.numeric(tapply(pSubs, .spacesPop, sum)), rep(1,max(.spacesPop)) )) )
			stop("divideTwDEMCSteps (1): pSubs within subspaces do not sum to 1.")
		#
		#.spacesPopOrig <- getSpacesPop(mcApp)
		#.pSpacesOrig <- as.numeric(tapply(pSubs, .spacesPopOrig, sum)
		boPopsDistConverge <- all(resStep$subPercChange <= subPercChangeCrit)
		mcApp0 <- squeeze(mcApp, length.out=ceiling(getNSamples(mcApp)*pThinPast) ) # mcApp0: thinned former mcApp
		mcApp <- mcApp1 <- .concatTwDEMCRuns(mcApp0,mc,doRecordProposals=doRecordProposals) # mcApp and mcApp1: thinned past + new samples 
		.nS1 <- sum(getNSamples(mcApp1))
		#
		if( nGenStep/thin >= 64){  # only check and split populations, if enough samples in batch
			#iPop=length(mc$pops)
			#---- check the populations for problems
			# TODO: think about plitting bevore sampling the first batch
			if( iGen != 0 ) for( iPop in seq_along(mc$pops)) if( nSamplePop[iPop] > 5){
						# see .tmp.f below for commands to trace error
						pop <- mc$pops[[iPop]] 
						#mtrace(.checkProblemsSpectralPop)
						resFCheck <- .checkProblemsSpectralPop( pop, critSpecVarRatio= critSpecVarRatio)
						if( resFCheck$hasProblems ){
							if( !is.null(dumpfileBasename) )
								if( dumpfileBasename == "recover"){
									# see .debug.divideTwDEMC below for debuggin commands
									cat("divideTwDEMCSteps: encountered convergence problems in a subSpace. calling recover() \n ")
									recover()
								}  else dump.frames(dumpfileBasename,TRUE)
							stop(paste("divideTwDEMCSteps: checking function encountered problems. Dump in ",dumpfileBasename,sep=""))
						}
					}
			#		
			#-------- split subspaces that are large and have changing quantiles
			iPopsSplit <- which( (resStep$pSubs == 1) | ((resStep$pSubs/2 > minPSub) & (resStep$subPercChange > subPercChangeCrit)) )
			if( 0 != length(iPopsSplit) ){
				#iPop = iCheckSplit[ length(iCheckSplit) ]
				newPops <- list()
				newPSubs <- integer(0)
				iPop <- iPopsSplit[1]
				for( iPop in iPopsSplit){
					popMc <- mc$pops[[iPop]]			# mc holds only the new samples
					aSample <- stackChains( popMc$parms )	# base splitting only on samples obtained in last batch
					#mtrace(getSubSpaces)
					subs <- getSubSpaces( aSample, isBreakEarly=FALSE, pSub=resStep$pSubs[iPop]
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
							stop("divideTwDEMCSteps: sample number does not match after splitting.")
						#.getParBoundsPops(c(newPopsI, list(pop)))
						#.getParBoundsPops(c(subs$spaces, list(pop))) #only internal splits
					}else{
						iPopsSplit[match(iPop,iPopsSplit)] <- NA		# remove from pops to split (do not remove in mcApp)
					}
				}
				iPopsSplit <- iPopsSplit[ is.finite(iPopsSplit) ] 
				if( 0 != length(iPopsSplit)){
					mcApp$pops <- c( mcApp$pops[-iPopsSplit], newPops )
					pSubs1 <- c( pSubs[-iPopsSplit], newPSubs)
				}
				# check consistency after splitting		
				lapply( mcApp$pops, .checkPop, nGen=12 )		
				.spacesPop <- getSpacesPop(mcApp)
				if( !isTRUE( all.equal(.pSpaces <- as.numeric(tapply(pSubs1, .spacesPop, sum)), rep(1,max(.spacesPop)))) )
					stop("divideTwDEMCSteps (2c): pSubs within subspaces do not sum to 1.")
				#
				.nS2 <- sum(getNSamples(mcApp))
				#if( (.nS2 > .nS1) || (.nS2 < .nS1-(nChainPop-1)*length(iPopsSplit)) )
				if( (.nS2 > .nS1) || (.nS2 < .nS1-(nChainPop-1)*(getNPops(mcApp)-getNPops(mcApp0))) )
					stop("divideTwDEMCSteps (2d): sample number does not match after splitting.")
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
							if( length(jPops)==0 ) stop("divideTwDEMCSteps: encountered no-neighboring case.")
							tmp.f <- function(){
								lapply(popsS, "[[", "splits")
								.spacePop <- getSpacesPop(mcApp0)
								lapply(subPops(mcApp0, iSpace=spaceInds[1])$pops, "[[", "splits")
								
							}
						}
				}
				#.getParBoundsPops(newPops)
			}
			#
			#-------- merge subspaces that contain less samples than minPSub/2
			# again base decision on sample of last batch (pSubs1)
			mcApp2 <- mcApp
			.nS2 <- sum(getNSamples(mcApp2))
			pSubs2 <- pSubs1
			mcApp <- mcApp2
			pSubs1 <- pSubs2
			iPopsMerge <- iPopsMerge0 <- which( pSubs1 < minPSub/2 | getNSamples(mcApp) < m0 )
			iExitMerge <- getNPops(mcApp2)	# prevent infinite runs
			.nMerges <- 0
			#if( 0 != length(iPopsMerge) ) recover()
			while( 0 != length(iPopsMerge) && iExitMerge != 0){
				#print(iPop)
				iPop <- iPopsMerge[ order(pSubs1[iPopsMerge])[1] ] # the one with lowest p
				#str3(mcApp$pops[[iPop]])
				#mcApp$pops[[iPop]]$splits
				#iPopsSpace <- which(getSpacesPop(mcApp) == mcApp$pops[[iPop]]$spaceInd)
				#lapply(mcApp$pops[iPopsSpace], "[[", "splits")
				#lapply(mcApp$pops[iPopsSpace], "[[", "splits")
				#mtrace(.mergePopTwDEMC)
				mcAppPrev <- mcApp; pSubsPrev <- pSubs1
				#mcApp <- mcAppPrev; pSubs1 <- pSubsPrev
				#mtrace(.mergePopTwDEMC)
				resMerge <- .mergePopTwDEMC( mcApp$pops, iPop, pSubs1 )
				# seek neighboring populations
				for( iSpace in 1:getNSpaces(mcApp) ){
					popsS <- resMerge$pops[ sapply(resMerge$pops, "[[","spaceInd") == spaceInds[iSpace] ]
					.nPopS <- length(popsS)
					.iPopsS <- seq_along(popsS)
					if( .nPopS > 1) for( iPop in  .iPopsS){
							.popSource <- popsS[[iPop]]
							jPops <- .iPopsS[-iPop][ sapply( popsS[-iPop], function(pop){ 
										(length(pop$splits) >= length(.popSource$splits)) && all(pop$splits[ 1:length(.popSource$splits)] == .popSource$splits)
									})]
							if( length(jPops)==0 ) stop("divideTwDEMCSteps_resMerge: encountered no-neighboring case.")
						}
				}
				#mcApp$pops[[iPop]]$spaceInd
				#rms <- sapply(resMerge$pops,"[[","spaceInd");	iPopsSpaceNew <- which(rms == 1); (tmp1<-lapply(resMerge$pops[iPopsSpaceNew], "[[", "splits"))
				#rms <- sapply(resMerge$pops,"[[","spaceInd");	iPopsSpaceNew <- which(rms == 2); (tmp2<-lapply(resMerge$pops[iPopsSpaceNew], "[[", "splits"))
				#print(sum(sapply(lapply(mcApp$pops,"[[","parms"),nrow)));	print(sum(sapply(lapply(resMerge$pops,"[[","parms"),nrow)))
				#
				#if( sum( sapply(lapply(resMerge$pops,"[[","parms"),nrow) ) > .nS2 ) { print("divideTwDEMCSteps: encountered larger sample size"); recover() }	
				mcApp$pops <- resMerge$pops
				pSubs1 <- resMerge$pPops
				iPopsMerge <- which( pSubs1 < minPSub/2 )
				tmp <- lapply( mcApp$pops, .checkPop, nGen=12 )		
				.spacesPop <- getSpacesPop(mcApp)
				if( !isTRUE( all.equal(.pSpaces <- as.numeric(tapply(pSubs1, .spacesPop, sum)), rep(1,max(.spacesPop)))) )
					stop("divideTwDEMCSteps (3): pSubs within subspaces do not sum to 1.")
				.nMerges <- .nMerges + resMerge$nSubs
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
							if( length(jPops)==0 ) stop("divideTwDEMCSteps_merge: encountered no-neighboring case.")
							tmp.f <- function(){
								lapply(popsS, "[[", "splits")
								.spacePop <- getSpacesPop(mcApp2)
								lapply(subPops(mcApp2, iSpace=spaceInds[1])$pops, "[[", "splits")
							}
						} #for iPop
				} #for iSpace
				#
				#print(iPopsMerge)
				iExitMerge <- iExitMerge -1 
			}
			if( iExitMerge == 0 ) stop("divideTwDEMCSteps: while loop of populations to merge did not exit.")
			# check for consistency
			.nS3 <- sum(sapply( mcApp$pops, function(jPop){ nrow(jPop$parms) }))
			if( (.nS3 > .nS2) || (.nS3 < .nS2-(nChainPop-1)*.nMerges) )
				stop("divideTwDEMCSteps (4): sample number does not match after merging.")
			if( any(getNSamples(mcApp) < m0) ){ print("Too few samples, recover"); recover() }
		} # enough new samples for checking/splitting
		
		iGen <- iGen + nGenStep
		#cat(paste(iGen," out of ",nGen," generations completed. T=",paste({T<-res$temp[nrow(res$temp),];round(T,digits=ifelse(T>20,0,1))},collapse="  "),"     ",date(),"\n",sep=""))
		cat(paste(iGen," out of ",nGen," generations completed. ,max(subPercChange)=",max(resStep$subPercChange),"     ",date(),"\n",sep=""))
	}	# while iGen < nGenBatch
	resStep$resTwDEMC <- mcApp
	resStep
}
attr(divideTwDEMCSteps,"ex") <- function(){
	data(den2dCorEx)
	#trace("twDEMCBlockInt", recover)
	mc0 <- den2dCorEx$mcSubspaces0
	#plot( as.mcmc.list(mc0), smooth=FALSE )
	#plot( 
	#mc0$pops[[ length(mc0$pops) ]]$TStart <- mc0$pops[[ length(mc0$pops) ]]$TEnd <- 100
	mc0$blocks[[1]]$fUpdateBlock <- updateBlockTwDEMC	# if beeing traced
	#mc0 <- res$resTwDEMC
	#mtrace(divideTwDEMCSteps)
	res <- divideTwDEMCSteps(mc0
		#, nGen=256*2+32
		, nGen=256*5
		, nGenBatch=256
		, dInfos=list(list(fLogDen=den2dCor))
		,  debugSequential=TRUE
		,  controlTwDEMC=list(DRgamma=0.1)
		, minPSub=0.05
	)
	getNSamples(res$resTwDEMC)
	sort(res$subPercChange, decreasing=TRUE)
	#str(controlTwDEMC)
	#str(list(...)$controlTwDEMC)	
	lapply(res$resTwDEMC$pops,"[[","splits")
	#.getParBoundsPops(res$resTwDEMC$pops)
	
	#windows(record=TRUE)
	plot( as.mcmc.list(stackPopsInSpace(mc0)), smooth=FALSE)
	#mc1 <- stackPopsInSpace(res$resTwDEMC, mergeMethod="stack")
	mc1 <- stackPopsInSpace(res$resTwDEMC, mergeMethod="random")
	plot( as.mcmc.list(mc1), smooth=FALSE) # note how the distribution shifted
	#plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=4)), smooth=FALSE)
	#plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=30)), smooth=FALSE)
	
	ss <- stackChains(concatPops(mc1))			# the new sample
	ss0 <- stackChains(concatPops(mc0))
	.nr <- max( nrow(ss), nrow(ss0) )
	ss <- ss[1:.nr,]
	ss0 <- ss0[1:.nr,]
	plot( b ~ a, as.data.frame(ss0), xlim=c(-0.5,2), ylim=c(-20,40) )
	plot( b ~ a, as.data.frame(ss), xlim=c(-0.5,2), ylim=c(-20,40) ) # not that more samples are in the region of interest
	points( 0.8,0.8, col="red" )
	plot( b ~ a, as.data.frame(ss), xlim=c(-2,2), ylim=80*c(-20,40) ) # not that more samples are in the region of interest
	points( 0.8,0.8, col="red" )
	
	gridx <- a <- seq(-0.5,2,length.out=91)
	gridy <- seq(-20,+40,length.out=91)
	gridX <- expand.grid(gridx, gridy)
	luden <- apply( gridX, 1, den2dCor ) 
	mLuden <- matrix(luden,nrow=length(gridx))
	image( gridx, gridy,  matrix(exp(luden),nrow=length(gridx)), col = rev(heat.colors(100)), xlab="a", ylab="b" )
	points( b ~ a, as.data.frame(ss), xlim=c(-0.5,2), ylim=c(-20,40) ) # not that more samples are in the region of interest
	points( 0.8,0.8, col="blue" )
	
}

