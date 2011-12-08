# see also subspace.R

divideTwDEMCStep <- function(
	### run twDEMCBlock on subspaces
	aTwDEMC					##<< former run of twDEMCBlockInt
	,spacePop=seq_along(aTwDEMC$pops)	##<< integer vector (nPop): specifying the space replicated that the population belongs to
	,qPop=numeric(0)		##<< numeric vector (nPop): probability of samples in subspace within entire space 
	,nGen=512				##<< the number of generations for twDEMCBlock
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMCBlock}} containing entry thin
	, ...					##<< further arguments to \code{\link{twDEMCBlock}}
){
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
		lapply(1:nPop, function(iPop){ cumNPopSpace[iPop]+1-(nPopSpace[iPop]:1)})	
	}
	nSamplePop <- getNSamples(aTwDEMC)
	
	#--  calculating initial quantiles and number of genrations
	if( 0 == length(qPop) ){
		#iSpace=nSpace
		nSampl
		qPopSpace <- lapply( iSpaces, function(iSpace){
				nSampleIPop <- nSamplePop[iPopsSpace[[iSpace]] ]
				nSampleIPop / sum(nSampleIPop)
			} ) 
		qPop <- do.call(c,qPopSpace)
	}
	nGen0 <- pmax(m0*thin, ceiling(minNSamplesSub*thin/nChainPop), nGen*qPop)		# at minimum m0*thin generations to keep sufficient samples for extending the run
	nGen0Thin <- (ceiling(nGen0/thin)+2)*thin		# +2 to avoid extending runs with little shift 
	
	# setting up initial populations (each sub of overall populations becomes one population)
	subs0 <- vector("list",nPop)	
	#iSub <- nSub
	for( iSub in iPops){
		subSpacesI <- subSpacesFlat[[iSub]]
		s0 <- subSpacesI$sample[sample.int(nrow(subSpacesI$sample),nInit), ]
		.tmp <- matrix( 1:(m0*nChainPop), ncol=nChainPop )
		#iChainPop <- 4
		Zinit <- abind( lapply(1:nChainPop, function(iChainPop){s0[.tmp[,iChainPop],]}), rev.along=0)
		names(dimnames(Zinit)) <- c("steps","parms","chains")
		subs0[[iSub]]$parms <- Zinit
		subs0[[iSub]]$upperParBounds <- subSpacesI$upperParBounds
		subs0[[iSub]]$lowerParBounds <- subSpacesI$lowerParBounds
	}
	
	resTwDEMC <- resTwDEMC0 <- twDEMCBlockInt(pops=subs0, nGen=nGen, controlTwDEMC=controlTwDEMC, ...)
	##details<< 
	## A first estimate of the proportion of samples from different subspaces
	## are the initial percentiles qp. 
	## The proportion of the samples from different subspaces is estimated
	## by the proportions of integrated density of the subspaces.
	## These proportions are esimated by average density over all samples multiplied by an
	## estimate proportional to the volume: the initial quantiles.
	#iPop <- 1
	#iPop <- 2
	popLogMeanDensSubs <- lapply( 1:nPop, function(iPop){
			iSubs <- iPopsSpace[[iPop]]
			#iSub <- 2
			#iSub <- nSub
			lw <- sapply(iSubs,function(iSub){  # average the unnormalized densities
					ss <- stackChains(resTwDEMC$pops[[iSub]]$logDen)
					ssLogDen <- rowSums(ss)
					twLogSumExp(ssLogDen, shiftUpperBound=TRUE)-log(length(ssLogDen)) 
				})	
			lw - max(lw, na.rm=TRUE)			#divide by a constant on exp scale to avoid numerical errors
		})	
	logMeanDensSubs <- do.call(c, popLogMeanDensSubs )
	#barplot(logMeanDensSubs)
	
	# estimate proportion of subspaces in the limiting distribution
	# it will be used to sample from the subspaces 
	pSubs <- do.call(c, lapply( 1:nPop, function(iPop){
				p2u <- exp(popLogMeanDensSubs[[iPop]])*qPopSpace[[iPop]]		# two weights: quantile and density importance
				p2 <- p2u/sum(p2u)									# normalize to 1 within population 
			}))
	subPercChange <- pSubs/qPop	# ratio of estimated proportion in limiting distribution to proportion of initial proportion 
	nSamplesSubsReq <- round(nGen/thin * pSubs)
	# because of rounding small differences in sum of sample numbers may occur 
	# for small deviations, adjust the number of samples in the largest sample
	# add argument nGen to pops
	for( iPop in 1:nPop ){
		dGen <- nGenThin/thin - sum(nSamplesSubsReq[iPopsSpace[[iPop]] ])  
		if( dGen != 0 ){
			if( abs(dGen)>length(iPopsSpace[[iPop]]) ) stop("divideTwDEMC: problem in calculating subspace sample numbers")
			iSubMax <- iPopsSpace[[iPop]][ which.max( nSamplesSubsReq[iPopsSpace[[iPop]] ]) ]
			nSamplesSubsReq[iSubMax] <- nSamplesSubsReq[iSubMax] + dGen
		} 
	}
	
	# extend the previous runs
	nSamplesAdd <- pmax(0, nSamplesSubsReq-getNSamples(resTwDEMC0))
	resTwDEMC <- resTwDEMC1 <- if( max(nSamplesAdd) == 0) resTwDEMC0 else
			#mtrace(twDEMCBlock.twDEMCPops)
			#mtrace(twDEMCBlockInt)
			resTwDEMC <- twDEMCBlock( resTwDEMC0, nGen=nSamplesAdd*thin, controlTwDEMC=controlTwDEMC, ...  )
	#all( getNSamples(resTwDEMC) >= nSamplesSubsReq)
	
	# do a subsampling
	# iPop=2
	subSamplePop <- function(iPop){
		iSubs  <- iPopsSpace[[iPop]] 
		#nSamplesSubsReq[iSubs]
		#(iSub <- iSubs[length(iSubs)])
		ssImpList <-  lapply( iSubs, function(iSub){
				if( nSamplesSubsReq[iSub] == 0) c() else{
					aTwDEMCSub <- concatPops(subPops( resTwDEMC, iPops=iSub ))
					ss <- stackChains(aTwDEMCSub)
					ssp <- ss[,-(1:getNBlocks(aTwDEMCSub))]
					# see .tmp.f below for commands to trace error
					if( nSamplesSubsReq[iSub] == 1) 
						ss[ sample.int( nrow(ss),1),] 
					else {
						resFCheck <- do.call( fCheckProblems, c(list(aTwDEMCSub), argsFCheckProblems))
						if( iSub==21) recover()
						if( resFCheck$hasProblems ){
							if( !is.null(dumpfileBasename) )
								if( dumpfileBasename == "recover"){
									cat("encountered convergence problems in a subSpace. calling recover() \n ")
									recover()
								}  else dump.frames(dumpfileBasename,TRUE)
							stop(paste("divideTwDEMC: checking function encountered problems. Dump in ",dumpfileBasename,sep=""))
						}
						sst <- ss[ round(seq(1,nrow(ss),length.out=nSamplesSubsReq[iSub])),]
					}
				}
			})
		ssImp <- do.call( rbind, ssImpList )
	}
	#mtrace(tmp.f)
	resl <- lapply( 1:nPop, subSamplePop)	# end lapply iPop: giving a single combined sample for each population
	#res <- ssImpPops <- abind(resl, rev.along=0)
	
	##value<< 
	## For each population, a list with entries
	res <- lapply(1:nPop,function(iPop){
			iSubs <- iPopsSpace[[iPop]]
			list(
				sample=resl[[iPop]]					##<< numeric matrix: the sampled parameters (rows: cases, cols: 1: logDensity, others parameter dimensions
				,subPercChange=subPercChange[iSubs]	##<< numeric vector: estimated proportion in limiting distribution to proportion of initial proportion 
				#,upperParBounds=subSpacesFlat[iSubs] upperParBounds[iSubs]	##<< list of named numeric vectors: the upper parameter bounds for the subspaces
				#,lowerParBounds=lowerParBounds[iSubs]	##<< dito upperParBounds
				,iVars=subSpacesL[[iPop]]$iVars		##<< reordered index or parameter dimensions to check for split, see \code{\link{findSplit}}
				,jVarsVar=subSpacesL[[iPop]]$jVarsVar		##<< reordered index or parameter dimensions to check for split, see \code{\link{findSplit}}
			#,pSubs0=pSubs0[posPop]
			#,pSubs2=p2[posPop]
			)
		})
	##end<<
}
attr(divideTwDEMCStep,"ex") <- function(){
	data(den2dCorTwDEMC)
	aTwDEMC0 <- thin(den2dCorTwDEMCPops, start=300)
	tmp <- squeeze(aTwDEMC0, length.out= 256 %/% getNChainsPop(aTwDEMC0) ) # thin to 256 samples per space
	ss0 <- stackChainsPop(concatPops(tmp))
	plot( b ~ a, as.data.frame(ss0[,,1]), xlim=c(-0.5,2), ylim=c(-20,40) )
	
	nBlock <- attr(ss0,"nBlock")
	minPSub <- 0.1
	subSpaces <- lapply( 1:dim(ss0)[3], function(iPop){
			samplePop <- ss0[,,iPop]
			#mtrace(getSubSpaces)
			#mtrace(findSplit)
			getSubSpaces(samplePop[,-(1:nBlock)], minPSub=minPSub, isBreakEarly=FALSE, argsFSplit=list())	# here omit the logDensity column
		})
	aTwDEMC <- dividePops( aTwDEMC0, subSpaces )
	
	
	
	
	
	res <- divideTwDEMCStep(aTwDEMC0)
	#aSamplePure <- aSample[,-1,]
	#aSamplePurePop <- aSamplePure[,,1]
	#mtrace(divideTwDEMC)
	#res <- divideTwDEMC(aSample, nGen=128, dInfos=list(list(fLogDen=den2dCor)),  isBreakEarly=FALSE, debugSequential=TRUE )
	res <- divideTwDEMC(aSample, nGen=256, dInfos=list(list(fLogDen=den2dCor)),  isBreakEarly=FALSE, debugSequential=TRUE,  controlTwDEMC=list(DRGamma=0.1) )
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
