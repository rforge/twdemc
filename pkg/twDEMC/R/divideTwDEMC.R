# see also subspace.R

divideTwDEMCStep <- function(
	### run twDEMCBlock on subspaces
	aTwDEMC					##<< former run of twDEMCBlockInt
	,qPop=numeric(0)		##<< numeric vector (nPop): probability of samples in subspace within entire space 
	,nGen=512				##<< the number of generations for twDEMCBlock
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMCBlock}} containing entry thin
	, doRecordProposals=FALSE		##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y
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
	iSpaces <- 1:nSpace
	iPopsSpace <- {						# index of pops in flat version for each space
		cumNPopSpace <- cumsum(nPopSpace)
		lapply(1:nSpace, function(iSpace){ cumNPopSpace[iSpace]+1-(nPopSpace[iSpace]:1)})	
	}
	nSamplePop <- getNSamples(aTwDEMC)
	
	#----  calculating initial quantiles and number of genrations
	if( 0 == length(qPop) ){
		#iSpace=nSpace
		qPopSpace <- lapply( iSpaces, function(iSpace){
				nSampleIPop <- nSamplePop[iPopsSpace[[iSpace]] ]
				nSampleIPop / sum(nSampleIPop)
			} ) 
		qPop <- do.call(c,qPopSpace)
	}
	#nGen0 <- pmax(m0*thin, ceiling(minNSamplesSub*thin/nChainPop), nGen*qPop)		# at minimum m0*thin generations to keep sufficient samples for extending the run
	nGen0Pops <- nGen*qPop		
	nGen0PopsThin <- (ceiling(nGen0Pops/thin))*thin		 
	
	#----- initial run based on current quantiles
	#trace("twDEMCBlockInt")
	resTwDEMC <- resTwDEMC0 <- twDEMCBlock(aTwDEMC, nGen=nGen0PopsThin, extendRun=FALSE, controlTwDEMC=controlTwDEMC, doRecordProposals=doRecordProposals, ...)
	recover()
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
	pPops <- do.call(c, lapply( 1:nSpace, function(iSpace){
				p2u <- exp(popLogMeanDensSubs[[iSpace]])*qPopSpace[[iSpace]]		# two weights: quantile and density importance
				p2 <- p2u/sum(p2u)									# normalize to 1 within population 
			}))
	subPercChange <- pPops/qPop	# ratio of estimated proportion in limiting distribution to proportion of initial proportion 
	nSamplesSubsReq <- round(nGenThin/thin * pPops)
	# because of rounding small differences in sum of sample numbers may occur 
	# for small deviations, adjust the number of samples in the largest sample
	# add argument nGen to pops
	#iSpace=nSpace
	for( iSpace in 1:nSpace ){
		dGen <- nGenThin/thin - sum(nSamplesSubsReq[iPopsSpace[[iSpace]] ])  
		if( dGen != 0 ){
			if( abs(dGen)>length(iPopsSpace[[iSpace]]) ) stop("divideTwDEMC: problem in calculating subspace sample numbers")
			iSubMax <- iPopsSpace[[iSpace]][ which.max( nSamplesSubsReq[iPopsSpace[[iSpace]] ]) ]
			nSamplesSubsReq[iSubMax] <- nSamplesSubsReq[iSubMax] + dGen
		} 
	}
	
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
			resTwDEMC <- twDEMCBlock( resTwDEMC0, nGen=nSamplesAdd*thin, controlTwDEMC=controlTwDEMC, doRecordProposals=doRecordProposals, ...  )
		}
	#all( getNSamples(resTwDEMC) >= nSamplesSubsReq)
	#getNSamples(resTwDEMC1)
	
	#----- do a subsampling
	resTwDEMC <- resTwDEMC2 <- squeeze.twDEMCPops(resTwDEMC1, length.out=nSamplesSubsReq)
	#getNSamples(resTwDEMC2)

	#----- append to former run
	#resTwDEMC <- .concatTwDEMCRuns(aTwDEMC,resTwDEMC2,doRecordProposals=doRecordProposals)
	#getNSamples(resTwDEMC)
	
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
	data(den2dCorTwDEMC)
	aTwDEMC <- den2dCorTwDEMCSpaces
	getSpacesPop(aTwDEMC)
	getNSamples(aTwDEMC)

	#mtrace(divideTwDEMCStep)
	res <- divideTwDEMCStep(aTwDEMC, nGen=256, dInfos=list(list(fLogDen=den2dCor))
		,  debugSequential=TRUE
		#,  controlTwDEMC=list(DRGamma=0.05) 
	)
	getNSamples(res$resTwDEMC)
	
	#windows(record=TRUE)
	plot( as.mcmc.list(stackPops(aTwDEMC)), smooth=FALSE)
	tmp <- stackPops(res$resTwDEMC, mergeMethod="stack")
	plot( as.mcmc.list(tmp), smooth=FALSE)
	plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=4)), smooth=FALSE)
	plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=30)), smooth=FALSE)
	plot( as.mcmc.list(subChains(concatPops(res$resTwDEMC,minPopLength=30),iPops=2)), smooth=FALSE)
	
	ss <- stackChains()
	plot( b ~ a, as.data.frame(res[[1]]$sample), xlim=c(-0.5,2), ylim=c(-20,40) )
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



divideTwDEMCSteps <- function(
	### run twDEMCBlock on subspaces
	aTwDEMC					##<< former run of twDEMCBlockInt
	,nGen=512				##<< the number of generations for twDEMCBlock
	, ...					##<< further arguments to \code{\link{twDEMCBlock}}
	,nGenBatch=256
	, argsFSplitPop=vector("list",dim(aSample)[3])	##<< for each population: list of arguments  passed \code{\link{getSubSpaces}} and further to \code{\link{findSplit}}, e.g. for passing order of variables to check in \code{iVars} and \code{jVarsVar}
	, dumpfileBasename="recover"
){
	nPop <- getNPops(aTwDEMC)
	resStep <- divideTwDEMCStep(aTwDEMC, nGen=nGenBatch, ... )
	
	#---- check the populations for problems
	mc <- resStep$resTwDEMC
	nSamplePop <- getNSamples(mc)
	iPop=nPop
	#for( iPop in 1:nPop) if( nSamplePop[iPop] > 1){
		mcp <- concatPops( mcpBlock <- subPops(mc, iPops=iPop ))
		ss <- stackChains(mcp)
		ssp <- ss[,-(1:getNBlocks(mcp))]
		
		recover()
		
		# see .tmp.f below for commands to trace error
		resFCheck <- do.call( fCheckProblems, c(list(mcp), argsFCheckProblems))
		#if( iPop==21) recover()
		if( resFCheck$hasProblems ){
			if( !is.null(dumpfileBasename) )
				if( dumpfileBasename == "recover"){
					cat("encountered convergence problems in a subSpace. calling recover() \n ")
					recover()
				}  else dump.frames(dumpfileBasename,TRUE)
			stop(paste("divideTwDEMC: checking function encountered problems. Dump in ",dumpfileBasename,sep=""))
		}
	#}

	##XXTODO: split pops with large proportions

	##XXTODO: merge pops with small proportions
	resStep
}
attr(divideTwDEMCSteps,"ex") <- function(){
	data(den2dCorTwDEMC)
	#mtrace(divideTwDEMCSteps)
	#trace("twDEMCBlockInt", recover)
	mc0 <- den2dCorTwDEMCSpaces
	mc0$pops[[ length(mc0$pops) ]]$T0 <- mc0$pops[[ length(mc0$pops) ]]$TEnd <- 100
	mc0$blocks[[1]]$fUpdateBlock <- updateBlockTwDEMC	# if beeing traced
	res <- divideTwDEMCSteps(mc0
		, nGen=256
		, nGenBatch=256
		, dInfos=list(list(fLogDen=den2dCor))
		,  debugSequential=TRUE
		,  controlTwDEMC=list(DRGamma=0.1)
	)
	getNSamples(res$resTwDEMC)
	#str(controlTwDEMC)
	#str(list(...)$controlTwDEMC)	
	
}

.debug.divideTwDEMC <- function(){
	str(mcp)
	plot(as.mcmc.list(mcp), smooth=FALSE)
	matplot( mcp$pAccept[,1,], type="l" )
	matplot( mcp$logDen[,1,], type="l" )
	matplot( mcp$temp, type="l" )
	plot(ss[,"a"], ss[,"b"], col=rev(heat.colors(100))[ twRescale(ss[,1],c(20,100)) ] )
	plot(mcp$Y[,"a",1], mcp$Y[,"b",1], col=rev(heat.colors(100))[ twRescale(mcp$Y[,"logDen1",1],c(20,100)) ] )
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


	





