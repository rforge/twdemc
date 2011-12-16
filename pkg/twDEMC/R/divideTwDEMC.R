# see also subspace.R

divideTwDEMCStep <- function(
	### run twDEMCBlock on subspaces
	aTwDEMC					##<< former run of twDEMCBlockInt
	,qPop=numeric(0)		##<< numeric vector (nPop): probability of samples in subspace within entire space 
	,nGen=512				##<< the number of generations for twDEMCBlock
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMCBlock}} containing entry thin
	, nSampleMin=16			##<< minimum number of samples in each population so that calculation of average density is stable
	, ...					##<< further arguments to \code{\link{twDEMCBlock}}
){
	spacePop=getSpacesPop(aTwDEMC)	# integer vector (nPop): specifying the space replicated that the population belongs to
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
	thin <- controlTwDEMC$thin
	nGenThin <- (nGen %/% thin)*thin		# when dealing with thinned samples use the number of generations of last recorded thinned generation
	nBlock <- getNBlocks(aTwDEMC)		# number of blocks in aTwDEMC
	nPop <- getNPops(aTwDEMC)			# number of populations
	iPops <- 1:nPop						
	nPopSpace <- table(spacePop)		# number of populaions per space replicated
	nSpace <- length(nPopSpace)			# number of space replicates
	spacesPop <- getSpacesPop(aTwDEMC)
	iPopsSpace <- lapply(1:nSpace, function(iSpace){ which(spacesPop==iSpace)}) # pops in flat version per space 
	nSamplePop <- getNSamples(aTwDEMC)
	
	#----  calculating initial quantiles and number of genrations
	if( 0 == length(qPop) ){
		#iSpace=nSpace
		qPopSpace <- lapply( 1:nSpace, function(iSpace){
				nSampleIPop <- nSamplePop[iPopsSpace[[iSpace]] ]
				nSampleIPop / sum(nSampleIPop)
			} ) 
		qPop <- do.call(c,qPopSpace)
	}
	#nGen0 <- pmax(m0*thin, ceiling(minNSamplesSub*thin/nChainPop), nGen*qPop)		# at minimum m0*thin generations to keep sufficient samples for extending the run
	nGen0Pops <- nGen*qPop
	m0 <- calcM0twDEMC(getNParms(aTwDEMC),getNChainsPop(aTwDEMC))		# minimum number of samples in step for extending runs
	nGen0PopsThin <- pmax( m0, nSampleMin, ceiling(nGen0Pops/thin))*thin	# must have at least m0 samples in order to be able to extend the sample	 
	
	#----- initial run based on current quantiles
	#trace("twDEMCBlockInt")
	resTwDEMC <- resTwDEMC0 <- twDEMCBlock(aTwDEMC, nGen=nGen0PopsThin, extendRun=FALSE, controlTwDEMC=controlTwDEMC, ...)
	#recover()
	#plot(as.mcmc.list(resTwDEMC),smooth=FALSE)
	#mcp <- concatPops( mcpBlock <- subPops(resTwDEMC0, iPops=nPop))
	#plot( as.mcmc.list(mcp), smooth=FALSE)
	
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
	popLogMeanDensSubs <- lapply( 1:nSpace, function(iSpace){
			iPops <- iPopsSpace[[iSpace]]
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
	
	# estimate proportion of subspaces in the limiting distribution
	# it will be used to sample from the subspaces 
	# iSpace=1
	pPops <- numeric( length(spacesPop))
	for( iSpace in 1:nSpace){
		p2u <- exp(popLogMeanDensSubs[[iSpace]])*qPopSpace[[iSpace]]		# two weights: quantile and density importance
		pPops[ iPopsSpace[[iSpace]] ] <- p2 <- p2u/sum(p2u)									# normalize to 1 within population 
	}
	subPercChange <- pPops/qPop	# ratio of estimated proportion in limiting distribution to proportion of initial proportion
	nSamplesSubsReq0 <- nGenThin/thin * pPops
	# because we sampled at least nSampleMin, we can append a factor of required samples
	
	# because of rounding small differences in sum of sample numbers may occur 
	# for small deviations, adjust the number of samples in the largest sample
	# add argument nGen to pops
	#iSpace=nSpace
	nSamplesSubsReq <- rep( NA_integer_, length(nSamplesSubsReq0) )
	for( iSpace in 1:nSpace ){ 
		extFac <- max(1,floor(min( (nGen0PopsThin/nSamplesSubsReq0)[ iPopsSpace[[iSpace]] ])))
		nSamplesSubsReq[ iPopsSpace[[iSpace]] ] <- nsReq <- round(extFac*nGenThin/thin * pPops[ iPopsSpace[[iSpace]] ])	# at least one sample in the subspace
	}
	#dGen <- extFac*nGenThin/thin - sum(nsReq)
	#	if( dGen != 0 ){
	#			if( abs(dGen)>length(iPopsSpace[[iSpace]]) ) stop("divideTwDEMC: problem in calculating subspace sample numbers")
	#		iSubMax <- iPopsSpace[[iSpace]][ which.max( nSamplesSubsReq[iPopsSpace[[iSpace]] ]) ]
	#		nSamplesSubsReq[iSubMax] <- nSamplesSubsReq[iSubMax] + dGen
	#	} 
	#}
	#nSamplesSubsReq <- pmax(1,nSamplesSubsReq)	# keep at least 1 sample per subspace
	
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
			resTwDEMC <- twDEMCBlock( resTwDEMC0, nGen=nSamplesAdd*thin, controlTwDEMC=controlTwDEMC,  ...  )
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
	pop			
	, critSpecVarRatio=50
){
	specVarRatio <- apply( pop$parms, 3, function(parmsChain){
			spec <- spectrum0.ar(parmsChain)$spec
			varSample <- apply(parmsChain, 2, var)
			spec/varSample
		})
	list(
		hasProblems = any(specVarRatio >= critSpecVarRatio)
		,specVarRatio = specVarRatio 
	)
}



divideTwDEMCSteps <- function(
	### run twDEMCBlock on subspaces
	aTwDEMC					##<< former run of twDEMCBlockInt
	,nGen=512				##<< the number of generations for twDEMCBlock
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMCBlock}} containing entry thin
	, doRecordProposals=FALSE		##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y
	, ...					##<< further arguments to \code{\link{twDEMCBlock}}
	,nGenBatch=256
	#, argsFSplitPop=vector("list",dim(aSample)[3])	##<< for each population: list of arguments  passed \code{\link{getSubSpaces}} and further to \code{\link{findSplit}}, e.g. for passing order of variables to check in \code{iVars} and \code{jVarsVar}
	, dumpfileBasename="recover"
	, critSpecVarRatio=20
	, minPSub = 0.1
	, subPercChangeCrit=1.6	##<< if all subPercChange of all sub-populations are below this value in two last batches, may assume convergence and skip further batches
	,maxNSample=128					##<< if given a value, then looking for subspaces is done on a subsample of given length for efficiency (see \code{\link{getSubSpaces}}) 
){
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
	thin <- controlTwDEMC$thin
	if( nGenBatch%/%thin < 64) 
		stop("not enough sample per batch. Increase nGenBatch to at least 64*thin or decreaste ctrolTwDEMC$thin")
	nGenBatch = (nGenBatch%/%thin)*thin		# make nGen multiple of thin
	boPopsDidConverge <- FALSE
	mcApp <- aTwDEMC		#new samples are appended each batch
	iGen = 0
	while( (iGen < nGen) && !boPopsDidConverge){
		
		#------- sample the next batch
		nGenStep <- min(nGenBatch, nGen-iGen)
		#mtrace(divideTwDEMCStep)
		resStep <- divideTwDEMCStep(mcApp, nGen=nGenStep, doRecordProposals=doRecordProposals, ... )
		mc <- resStep$resTwDEMC
		pSubs <- pSubs1 <- resStep$pSubs	#pSubs1 will reflect splitted populations
		#
		nSamplePop <- getNSamples(mc)
		.tmp <- lapply( mc$pops[nSamplePop != 0], .checkPop, nGen=12 )		
		.spacesPop <- getSpacesPop(mc)
		#tapply(nSamplePop, .spacesPop, function(nsI){ nsI/sum(nsI)})
		if( !isTRUE( all.equal(.pSpaces <- as.numeric(tapply(pSubs, .spacesPop, sum)), rep(1,max(.spacesPop)) )) )
			stop("divideTwDEMCSteps (1): pSubs within subspaces do not sum to 1.")
		
		#.spacesPopOrig <- getSpacesPop(mcApp)
		#.pSpacesOrig <- as.numeric(tapply(pSubs, .spacesPopOrig, sum)
		boPopsDistConverge <- all(resStep$subPercChange <= subPercChangeCrit)
		# XXTODO: think about thinning mcApp before appending new samples to downweigh history
		mcApp <- .concatTwDEMCRuns(mcApp,mc,doRecordProposals=doRecordProposals)
		lapply( mcApp$pops, .checkPop, nGen=12 )		
		
		
		if( nGenStep/thin >= 64){  # only check, split and merge populations, if enough samples in batch
			#iPop=length(mc$pops)
			#---- check the populations for problems
			for( iPop in seq_along(mc$pops)) if( nSamplePop[iPop] > 5){
				#mcp <- concatPops(  subPops(mc, iPops=iPop ))
				#plot( as.mcmc.list(mcp), smooth=FALSE )
				#matplot( mcp$pAccept[,1,], type="l" )
				#.getParBoundsPops(list(pop))
				#ss <- stackChains(mcp)
				#ssp <- ss[,-(1:getNBlocks(mcp))]
				# see .tmp.f below for commands to trace error
				pop <- mc$pops[[iPop]] 
				resFCheck <- .checkProblemsSpectralPop( pop, critSpecVarRatio= critSpecVarRatio)
				#if( iPop==21) recover()
				if( resFCheck$hasProblems ){
					if( !is.null(dumpfileBasename) )
						if( dumpfileBasename == "recover"){
							cat("encountered convergence problems in a subSpace. calling recover() \n ")
							recover()
						}  else dump.frames(dumpfileBasename,TRUE)
					stop(paste("divideTwDEMC: checking function encountered problems. Dump in ",dumpfileBasename,sep=""))
				}
			}
		
			#-------- split subspaces that are large and have changing quantiles
			iPopsSplit <- which( (resStep$pSubs/2 > minPSub) & (resStep$subPercChange >  subPercChangeCrit))
			if( 0 != length(iPopsSplit) ){
				#iPop = iCheckSplit[ length(iCheckSplit) ]
				newPops <- list()
				newPSubs <- integer(0)
				for( iPop in iPopsSplit){
					popMc <- mc$pops[[iPop]]
					aSample <- stackChains( popMc$parms )	# base splitting only on samples obtained in last batch
					#mtrace(getSubSpaces)
					subs <- getSubSpaces( aSample, isBreakEarly=FALSE, pSub=resStep$pSubs[iPop], minPSub=minPSub, maxNSample=maxNSample
					, splits=popMc$splits, lowerParBounds=popMc$lowerParBounds, upperParBounds=popMc$upperParBounds )
					#.getParBoundsPops(c(list(popMc),subs$spaces))
					if( length(subs$spaces) > 1){
						pop <- mcApp$pops[[iPop]]
						#mtrace(divideTwDEMCPop)
						newPopsI <- divideTwDEMCPop(pop, subs$spaces)
						newPops <- c( newPops, newPopsI)
						newPSubs <- c( newPSubs, sapply(subs$spaces, "[[", "pSub") )
						lapply( newPopsI, .checkPop, nGen=12 )
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
				#.getParBoundsPops(newPops)
			}
			
			#-------- merge subspaces that contain less samples than minPSub/2
			# again base decision on sample of last batch (pSubs1)
			iPopsMerge <- which( pSubs1 < minPSub/2 )
			iExitMerge <- getNPops(mcApp)	# prevent infinite runs
			while( 0 != length(iPopsMerge) && iExitMerge != 0){
				iPop <- iPopsMerge[ order(pSubs1[iPopsMerge])[1] ] # the one with lowest p
				resMerge <- .mergePopTwDEMC( mcApp$pops, iPop, pSubs1 )
				mcApp$pops <- resMerge$pops
				pSubs1 <- resMerge$pPops
				iPopsMerge <- which( pSubs1 < minPSub/2 )
				
				lapply( mcApp$pops, .checkPop, nGen=12 )		
				.spacesPop <- getSpacesPop(mcApp)
				if( !isTRUE( all.equal(.pSpaces <- as.numeric(tapply(pSubs1, .spacesPop, sum)), rep(1,max(.spacesPop)))) )
					stop("divideTwDEMCSteps (3): pSubs within subspaces do not sum to 1.")
				iExitMerge <- iExitMerge -1 
			}
			if( iExitMerge == 0 ) stop("divideTwDEMCSteps: while loop of populations to merge did not exit.")
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
	#mc0$pops[[ length(mc0$pops) ]]$T0 <- mc0$pops[[ length(mc0$pops) ]]$TEnd <- 100
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
	ss0 <- stackChains(concatPops(mc0))[1:nrow(ss0),]
	plot( b ~ a, as.data.frame(ss0), xlim=c(-0.5,2), ylim=c(-20,40) )
	plot( b ~ a, as.data.frame(ss), xlim=c(-0.5,2), ylim=c(-20,40) ) # not that more samples are in the region of interest
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

.debug.divideTwDEMC <- function(){
	str(mcp)
	plot(as.mcmc.list(mcp), smooth=FALSE)
	matplot( mcp$pAccept[,1,], type="l" )
	matplot( mcp$logDen[,1,], type="l" )
	matplot( mcp$temp, type="l" )
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


	





