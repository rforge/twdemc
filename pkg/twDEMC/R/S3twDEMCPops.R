#----------------------------------- subset ---------------
setMethodS3("subset","twDEMCPops", function( 
		### Condenses an twDEMCPops result object to the cases boKeep.
		x 			##<< twDEMCPops object
		,boKeep		##<< either logical vector or numeric vector of indices of cases to keep
		,...
		,iPops=seq_along(x$pops)	##<< integer vector: only these populations are subset, others are kept
		,dropShortPops=FALSE		##<< if set to TRUE, pops in iPops with less samples than to what \code{boKeep} refers to are dropped
	){
		# subset.twDEMCPops 
		
		##seealso<<   
		## \code{\link{twDEMCBlockInt}}
		
		##details<< 
		## There are several methods access properties a result of an \code{\link{twDEMCBlockInt}}, i.e.
		## an object of class \code{twDEMCPops} \itemize{
		## \item{ number of generations: \code{\link{getNGen.twDEMCPops}}  } 
		## \item{ number of samples (only one sample each thinning inteval): \code{\link{getNSamples.twDEMCPops}}  } 
		## \item{ number of chains: \code{\link{getNChains.twDEMCPops}}  } 
		## \item{ number of populations: \code{\link{getNPops.twDEMCPops}}  } 
		## \item{ number of chains per population: \code{\link{getNChainsPop.twDEMCPops}}  } 
		## \item{ number of parameters: \code{\link{getNParms.twDEMCPops}}  } 
		## \item{ thinning interval: \code{res$thin}  } 
		## \item{ space replicate that the poplation belongs to: \code{\link{getSpacesPop.twDEMCPops}}  } 
		## \item{ number of space replicates: \code{\link{getNSpaces.twDEMCPops}}  } 
		## \item{ number of blocks: \code{\link{getNBlocks.twDEMCPops}}  } 
		## \item{ parameter bounds: \code{\link{getParBoundsPop.twDEMCPops}}  } 
		## \item{ parameter bounds: \code{\link{getCurrentTemp.twDEMCPops}}  } 
		##}
		##
		## There are several methods to transform or subset the results of an \code{\link{twDEMCBlockInt}} run. \itemize{
		## \item{ transforming to array representation across populations of type \code{twDEMC}: code{\link{concatPops.twDEMCPops}}  }
		## \item{ select sub-populations: \code{\link{subPops.twDEMCPops}}  } 
		## \item{ select subset of cases: \code{\link{subset.twDEMCPops}} (this function) }
		## \item{ make populations the same length: \code{\link{squeeze.twDEMCPops}}  }
		## \item{ stack all the results of all chains to one big matrix: \code{\link{stackChains.twDEMCPops}}  } 
		## \item{ thin all the chains: \code{\link{thin.twDEMCPops}}  } 
		## \item{ combine several twDEMCPops results to a bigger set of populations: \code{\link{combineTwDEMCPops}}  }
		## \item{ combine populations for subspaces to bigger populations: \code{\link{stackPopsInSpace.twDEMCPops}}  }
		## \item{ combine MarkovChains of each population of a twDEMC to a single chain: \code{\link{stackChainsPop.twDEMCPops}}  }
		##}
		##
		## There are several methods utilize the functions of the coda package. \itemize{
		## \item{ convert an twDEMCPops to a coda mcmc.list \code{\link{as.mcmc.list.twDEMCPops}}  } 
		## \item{ applying a function to all of the chains: \code{\link{mcmcListApply}}  }
		## \item{ stack all the results of all chains to one big matrix: \code{\link{stackChains.mcmc.list}}  } 
		## \item{ transforming parameters \code{\link{transOrigPopt.mcmc.list}}  } 
		##}
		#? \item{ stack all the results of all chains to one big matrix \code{\link{stackChains.twDEMCPops}}  } 
		#? \item{ plotting a subset of the chains and cases: \code{\link{plotThinned.mcmc.list}}  } 
		
		
		# To help re-initializing the arguments to fLogDen \itemize{
		# \item{ transforming parameters \code{\link{subsetArgsFLogDen}}  }}
		
		
		nSamplesPop <- getNSamples(x)
		iKeep <- if( is.numeric(boKeep)) boKeep else which(boKeep)
		maxStep <- max(iKeep)
		if( dropShortPops ){
			x <- subPops( x, which( (1:seq_along(x$pops) == iPops) && (nSamplesPop >= maxStep)) )
		}else if( maxStep > min(nSamplesPop[iPops])) stop(
				"subset.twDEMCPops: provided indices outside the steps of the smalles population. Use dropShortPops=TRUE to drop shorter populations before.")
		for( iPop in iPops ){
			#mtrace(.subsetTwDEMCPop)
			x$pops[[iPop]] <- tmp <- .subsetTwDEMCPop(x$pops[[iPop]], iKeep)
		}
		##details<<
		## components \code{thin,Y,nGenBurnin} are kept, but may be meaningless after subsetting.
		x
		### list of class twDEMCPops with subset of cases in parsm, logDen, pAccept, and temp
	})
#mtrace(subset.twDEMCPops)
attr(subset.twDEMCPops,"ex") <- function(){
	if( FALSE){
		tmp <- subset(res,1:3)
		str(tmp)
		tmp <- subset(res,1:10)  # should produce an error
		#mtrace(subset.twDEMCPops)
		tmp <- subset(res,1:10,dropShortPops=TRUE)  
	} 
}



setMethodS3("getNGen","twDEMCPops", function( 
		### Extract the number of completed generations in res
		x	##<< object of class twDEMCPops
		,... 
	){
		##details<< 
		## the number of generations corresponds to the steps after time 0 (sample 1).
		## Hence the sample of size 2 and thinning 1 describes one generation (one step forward).
		## A sample of size 2 of thinning 5 encompasses times 0 and 5, i.e. 5 generations.
		## see 
		
		# getNGen.twDEMCPops
		##seealso<<   
		## \code{\link{iSample2time}}
		## \code{\link{time2iSample}}
		## \code{\link{getNPops.twDEMCPops}}
		## \code{\link{getNSamples.twDEMCPops}}
		## \code{\link{getNSamplesSpace.twDEMCPops}}
		## \code{\link{getNChains.twDEMCPops}}
		## \code{\link{getNChainsPop.twDEMCPops}}
		## \code{\link{getNParms.twDEMCPops}}
		## \code{\link{getNSpaces.twDEMCPops}}
		## \code{\link{getSpacesPop.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		#mtrace(getNSamples.twDEMCPops)
		nS <- getNSamples(x)
		ifelse( nS == 0, 0, (nS-1)*x$thin )
		### integer vector, number of completed generations
	})
attr( getNGen.twDEMCPops, "ex") <- function(){
	data(twdemcEx1)
	getNGen(twdemcEx1)
	getNSamples(twdemcEx1)
	getNSamplesPop(twdemcEx1)
	getNSamplesSpace(twdemcEx1)
	twdemcEx1$thin
	getNPops(twdemcEx1)
	getNChains(twdemcEx1)
	getNChainsPop(twdemcEx1)
	getNParms(twdemcEx1)
	getNBlocks(twdemcEx1)
}

setMethodS3("getNPops","twDEMCPops", function( 
		### Extracts the number of populations
		x	##<< object of class twDEMCPops
		,... 
	){
		# getNPops.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## ,\code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		length(x$pops)
		### integer, number of populations in twDEMCPops
	})

setMethodS3("getNChains","twDEMCPops", function( 
		### Extracts the number of chains
		x	##<< object of class twDEMCPops
		,... 
	){
		# getNChains.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		length(x$pops) * dim(x$pops[[1]]$parms)[3]
		### integer, number of chains in twDEMCPops
	})

setMethodS3("getNChainsPop","twDEMCPops", function( 
		### Extracts the number of chains per population
		x	##<< object of class twDEMCPops
		,... 
	){
		# getNChainsPop.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		dim(x$pops[[1]]$parms)[3]
		### integer, number of chains per population in twDEMCPops
	})

setMethodS3("getNParms","twDEMCPops", function( 
		### Extracts the number of parameters, i.e. the dimension of the parameter vector
		x	##<< object of class twDEMCPops
		,... 
	){
		# getNParms.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		ncol(x$pops[[1]]$parms)
		### integer, number of parameters in twDEMCPops
	})

setMethodS3("getNSamples","twDEMCPops", function( 
		### Extracts the number of samples
		x	##<< object of class twDEMCPops
		,... 
	){
		# getNSamples.twDEMCPops
		##seealso<<   
		## \code{\link{getNSamplesSpace.twDEMCPops}}
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		##details<< There is only one sample per thinning interval of length \code{x$thin}.
		sapply( x$pops, function(pop){ nrow(pop$parms) })
		### integer vector: number of samples in each population of twDEMCPops
	})

setMethodS3("getNSamplesSpace","twDEMCPops", function( 
		### Extracts the number of samples summed over population within one space
		x	##<< object of class twDEMCPops
		,... 
	){
		# getNSamples.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{getNSamples.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		##details<< There is only one sample per thinning interval of length \code{x$thin}.
		.nSamplesPop <- getNSamples(x)
		.spacePop <- getSpacesPop(x)
		spaceInds <- unique(.spacePop)
		.nSamplesSpace <- structure( sapply(spaceInds, function(spaceInd){ sum(.nSamplesPop[.spacePop==spaceInd])}), names=spaceInds)
		.nSamplesSpace
		### integer vector: number of samples in each population of twDEMCPops
	})


setMethodS3("getSpacesPop","twDEMCPops", function( 
		### Extracts the space replicate that each population belongs to.
		x	##<< object of class twDEMCPops
		,... 
	){
		# getSpacesPop.twDEMCPops
		##seealso<<   
		## \code{\link{getNSpaces.twDEMCPops}}
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		sapply( x$pops, "[[" ,"spaceInd" )
		### integer vector (nPop): index of the replicated space.
	})

setMethodS3("getNSpaces","twDEMCPops", function( 
		### Extracts the number of space replicates.
		x	##<< object of class twDEMCPops
		,... 
	){
		# getNSpaces.twDEMCPops
		##seealso<<   
		## \code{\link{getSpacesPop.twDEMCPops}}
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		length(unique(getSpacesPop.twDEMCPops(x)))
		### scalar integer: the number of spaces
	})


setMethodS3("getNBlocks","twDEMCPops", function( 
		### Extracts the number of samples
		x	##<< object of class twDEMCPops
		,... 
	){
		# getNBlocks.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		##details<< There is only one sample per thinning interval of length \code{x$thin}.
		length(x$blocks)
		### integer vector: number of samples in each population of twDEMCPops
	})

setMethodS3("getCurrentTemp","twDEMCPops", function( 
		### Get the Temperature, i.e. cost reduction factor of the last sample
		x	##<< object of class twDEMCPops
		,... 
	){
		# getCurrentTemp.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		##details<< 
		iPop <- which.max(getNSamples(x))	# very short pops may have not sampled low temperature
		x$pops[[iPop]]$temp[ nrow(x$pops[[iPop]]$temp), ]
		### numeric vector: Temperature for each result component for the last sample
	})



.getParBoundsPops <- function(
	pops
){
#	lapply( c("upperParBounds","lowerParBounds","splits"), function(entry){
#			sapply(pops, function(pop){ pop[entry] })
#		} )
	lapply(pops, function(pop){ pop[c("upperParBounds","lowerParBounds","splits")] })
}

setMethodS3("getParBoundsPop","twDEMCPops", function( 
		### Extracts lower or upper ParBounds
		x	##<< object of class twDEMCPops
		,... 
	){
		# getParBoundsPop.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCBlockInt}}
		
		##value<<
		list( ##describe<<
			upperParBoundsPop = lapply(x$pops, "[[", "upperParBounds")	##<< list of named numeric vectors giving upper parameter bounds per population 
			,lowerParBoundsPop = lapply(x$pops, "[[", "lowerParBounds") ##<< list of named numeric vectors giving lower parameter bounds per population
		) ##end <<
	})



#------------------------ concatPops -------------------------------------
setMethodS3("concatPops","twDEMCPops", function( 
	### Concatenates all the chains of all subpopulations to one matrix of call \code{twDEMC}.
	x
	,... 
	, useThinning=TRUE	##<< if TRUE thinning is used to make populations the same length, if FALSE they are cut to shortest population
	, minPopLength=NULL	##<< integer scalar: if specified, populations with less samples than length.out are dropped
){
	#concatPops.twDEMCPops
	##seealso<<   
	## \code{\link{subset.twDEMCPops}}
	## ,\code{\link{subset.twDEMC}}
	## ,\code{\link{twDEMCBlockInt}}
	nStepsPop <- getNSamples(x)
	if( 1 == length(minPopLength) ){
		iKeepPops <- which( nStepsPop >= minPopLength )
		x <- subPops(x, iPops=iKeepPops )
		nStepsPop <- nStepsPop[iKeepPops]
	}
	nSteps <- min(nStepsPop)
	if( nSteps < 2 ) stop("concatPops.twDEMCPops: min(nStepsPop)<2: cannot reduce to single state. use minPopLength=2 to drop degenerated populations")
	if( !all(nStepsPop == nSteps) ){
		if( useThinning)
			x <- squeeze(x, length.out=nSteps )
		else
			x <- subset(x, 1:nSteps)
	}
	pops <- x$pops
	p1 <- pops[[1]]
	x$pops <- NULL
	x$parms <- structure( abind( lapply(pops,"[[","parms"), along=3), dimnames=dimnames(p1$parms))
	#x$temp <- structure( abind( lapply(pops,"[[","temp"), along=3), dimnames=dimnames(p1$temp))
	x$temp <- p1$temp	# temperature of all population and chains have equal start and end
	x$pAccept <- structure( abind( lapply(pops,"[[","pAccept"), along=3), dimnames=dimnames(p1$pAccept))
	x$resLogDen <- structure( abind( lapply(pops,"[[","resLogDen"), along=3), dimnames=dimnames(p1$resLogDen))
	x$logDen <- structure( abind( lapply(pops,"[[","logDen"), along=3), dimnames=dimnames(p1$logDen))
	# for YL constrain to the 
	YL <-  lapply(pops,"[[","Y")
	nY <- min(sapply(YL,nrow))
	YLs <- lapply(YL, function(Y){ Y[nrow(Y)-((nY-1):0),,,drop=FALSE] })
	x$Y <- structure( abind(YLs, along=3), dimnames=dimnames(p1$Y))
	x$upperParBoundsPop = lapply( pops, "[[", "upperParBounds" )
	x$lowerParBoundsPop = lapply( pops, "[[", "lowerParBounds" )
	x$nPop <- length(pops)
	class(x) <- c("list","twDEMC")
	### An object of class twDEMC
	x
})
attr(concatPops,"ex") <- function(){
	if( FALSE ){
		getNSamples(tmp <- concatPops(res))
		getNChains(tmp)
		getNPops(tmp)
		#mtrace(concatPops.twDEMCPops)
		getNSamples(tmp <- concatPops(res,minPopLength=10))
		getNChains(tmp)
		getNPops(tmp)
	}
}

setMethodS3("subsetF","twDEMCPops", function( 
		### keeps only cases within population which evaluate to TRUE for a given function.
		x		##<< object of class twDEMCPops
		,fKeep	##<< function(pop) returning an boolean matrix (nStep x nChain) of cases to keep
			##<< ,alternatively returning an integer matrix (niStep x nChain) with the indices to keep
			##<< ,alternatively returning an integer or boolean vector, that is applied to each chain
		,...
	){
		#subsetF.twDEMCPops
		##details<< The samples are redistributed across chains
		xNew <- x
		nChainPop <- getNChainsPop(x)
		xNew$pops <- lapply(x$pops, function(pop){
				if( nrow(pop$parms) == 0){
					pop # no rows to apply fKeep to
				}else{		
					boKeep <- fKeep(pop)
					if( is.matrix(boKeep) ){
						# create a subset population for each chain
						.subSamplesChain <- lapply( 1:nChainPop, function(iChain){
								.subsetTwDEMCPop(pop, boKeep[,iChain], iChain)
							})
						# combine all the populations to a common population of one chain
						ss <- combineTwDEMCPops(.subSamplesChain, mergeMethod="random")			
						# split into pops, of equal length
						ssu <- .unstackPopsTwDEMCPops(ss$pop, nChainPop)	# here may loose or gain a few samples due to not a multiple of nChains
						# combine the pops of the chains into one populatin again
						newPop <- .concatChainsTwDEMCPops(ssu)	# combine chains into one array
						newPop$splits <- pop$splits
						newPop
					}else{
						# if only a vector all chains have the same length 
						# not need of combine/unstack
						.subsetTwDEMCPop(pop, boKeep)
					}# if is.matix 
				}# nrow(parms) != 0
			})
		### twDEMCPops with each population with some cases removed.
		xNew
	})
attr(subsetF.twDEMCPops,"ex") <- function(){
	data(twdemcEx1)
	range(concatPops(twdemcEx1)$parms[,"a",]) # spanning 9 to 11
	#pop <- twdemcEx1$pops[[1]]
	# note that the numer of samples across chains within one population is allowed to differ
	fKeep <- function(pop){ tmp <- (pop$parms[,"a",] < 10) }
	res <- subsetF(twdemcEx1, fKeep )
	plot( as.mcmc.list(res), smooth=FALSE )
	getNSamples(res)
}

setMethodS3("subsetTail","twDEMCPops", function( 
		### discards the first part of all the chains
		x		##<< object of class twDEMCPops
		,pKeep=0.5	##<< the percentage of the samples to keep 
		,... 
	){
		fKeep <- function(pop){
			nR <- nrow(pop$parms)
			nDrop <- round(nR*(1-pKeep))
			if( nDrop >= nR) return(FALSE) else return( (nDrop+1):nR )
		}
		subsetF(x,fKeep)
	})
attr(subsetTail.twDEMCPops,"ex") <- function(){
	data(twdemcEx1)
	res <- subsetTail(twdemcEx1)
	plot( as.mcmc.list(res), smooth=FALSE )
	getNSamples(twdemcEx1)
	getNSamples(res)
	
	# special case: no rows in one pop of twdemcEx1
	twdemcEx1$pops[[2]] <- .subsetTwDEMCPop
}




#mtrace(concatPops.twDEMCPops)
.subsetTwDEMCPop <- function(
	### subset all items of a twDEMCPops population
	pop		##<< single populatin of a twDEMCPops object
	,iKeep	##<< indices to keep (integer vector or boolean)
	,iChain=TRUE	##<< indices of chains to keep
){
	##details<< no checking of bounds performed
	newPop <- pop	# keep entries upperParBounds, lowerParBounds, splits
	if( nrow(pop$parms)==0 ) iKeep<-0
	newPop$parms <- pop$parms[iKeep,,iChain  ,drop=FALSE] 
	newPop$logDen <- pop$logDen[iKeep,,iChain ,drop=FALSE] 
	newPop$resLogDen <- pop$resLogDen[iKeep,,iChain ,drop=FALSE] 
	newPop$pAccept <- pop$pAccept[iKeep,,iChain ,drop=FALSE]
	newPop$temp <- pop$temp[iKeep, ,drop=FALSE]
	newPop
}

setMethodS3("subPops","twDEMCPops", function( 
		### Condenses an twDEMCPops List to given population indices or space replicates
		x
		, iPops 		##<< integer vector listing those spaceInd to keep
		, iSpaces=NULL	##<< alternatively specifying the space replicates for which to keep populations
		,... 
	#, doKeepBatchCall=FALSE	##<< wheter to retain the batch call attribute of x
	){
		# subPops.twDEMCPops
		##seealso<<   
		## \code{\link{twDEMCBlockInt}}
		## \code{\link{subset.twDEMCPops}}
		if( 0 != length(iSpaces) )
			iPops <- which(getSpacesPop(x) %in% iSpaces )
		x$pops <- x$pops[iPops]
		x
	})
#mtrace(subPops.twDEMCPops)



setMethodS3("squeeze","twDEMCPops", function(
	### Reduces the rows of populations so that all populations have the same number of samples. 
	x, ##<< the twDEMCPops list to thin 
	...,
	length.out=min(getNSamples(x)),	##<< number vector (nPops) of steps in each population, or numeric scalar specifying the lenght for all populations
	dropShortPops=FALSE				##<< if set to TRUE, pops with less samples than length.out[1] are dropped
){
	# squeeze.twDEMCPops
	##seealso<<   
	## \code{\link{subset.twDEMCPops}}
	nPop <- getNPops(x)
	nSamplesPop <- getNSamples(x)
	if( 1==length(length.out)){
		if( dropShortPops ){
			x <- subPops( x, iPops=which(nSamplesPop < length.out))
		}else if( length.out > min(nSamplesPop)) stop(
				"squeeze.twDEMCPops: specified a length that is longer than the shortest population. Use dropShortPops=TRUE to drop shorter populations.")
		length.out = rep( length.out, nPop)
	}
	if( nPop != length(length.out)) stop("squeeze.twDEMCPops: length.out must specify a length for each population.")
	length.out <- pmin( length.out, nSamplesPop ) # avoid lengthening the pops
	if( all( length.out == nSamplesPop) )
		return(x)
	##details<< 
	## all populations with 
	for( iPop in seq_along(x$pops) ){
		iKeep <- seq(1,nSamplesPop[iPop],length.out=length.out[iPop])
		iKeep[length(iKeep)] <- nSamplesPop[iPop]	# take the last row
		x$pops[[iPop]]$parms <- x$pops[[iPop]]$parms[iKeep,, ,drop=FALSE] 
		x$pops[[iPop]]$logDen <- x$pops[[iPop]]$logDen[iKeep,, ,drop=FALSE] 
		x$pops[[iPop]]$resLogDen <- x$pops[[iPop]]$resLogDen[iKeep,, ,drop=FALSE] 
		x$pops[[iPop]]$pAccept <- x$pops[[iPop]]$pAccept[iKeep,, ,drop=FALSE]
		x$pops[[iPop]]$temp <- x$pops[[iPop]]$temp[iKeep, ,drop=FALSE]
	}
	##details<<
	## components \code{thin,Y,nGenBurnin} are kept, but may be meaningless after subsetting.
	x
})
attr(squeeze.twDEMCPops,"ex") <- function(){
	if( FALSE){
		tmp <- squeeze(res)
		getNSamples(tmp)
		tmp2 <- subPops(res,2)
		#mtrace(squeeze.twDEMCPops)
		getNSamples( squeeze(tmp2,length.out=10) )
		getNSamples( squeeze(tmp,length.out=10) )	# should produce error
		getNSamples( squeeze(tmp,length.out=10, dropShortPops=TRUE) )	# should produce error
	} 
}

setMethodS3("stackChains","twDEMCPops", function( 
		### Combine logLik and parms of MarkovChains of a twDEMCPops object to one matrix. 
		x
		,...				##<< not used 
	){
		# stackChains.twDEMCPops
		##seealso<<   
		## \code{\link{stackPopsInSpace}}, \code{\link{subset.twDEMCPops}} 
		cPop <- combineTwDEMCPops(x$pops)$pop
		cArr <- abind( cPop$logDen, cPop$parms, along=2)
		res <- stackChains.array( cArr )
		### Matrix with first nDen columns the logDensity logDen and the remaining columns the variables.
	})

setMethodS3("thin","twDEMCPops", function( 
		### Reduces the rows of an twDEMCPops object (list returned by \code{\link{twDEMCBlockInt}}) to correspond to a thinning of \code{newThin}.
		x, ##<< the twDEMCPops list to thin 
		newThin=x$thin, ##<< finite numeric scalar: the target thinning factor, must be positive multiple of x$thin
		start=0,  
		### numeric scalar: the start time of the chain. 
		### Note that time starts from zero.
		### If a vector or matrix is supplied (e.g. nGenBurnin) then the maximum is used
		end=NA,   
		### numeric vector the maximum end time of the populations. 
		### Note that time starts from zero.
		### If a scalar is supplied, then it is repeateed NPop times 
		...
		, doKeepBatchCall=FALSE	##<< wheter to retain the batch call attribute of x
	
	){
		# thin.twDEMCPops
		##seealso<<   
		## \code{\link{subset.twDEMCPops}}
		
		# with the thinned list having mZ rows, this corresponds to (mZ-1)*thin Metropolis steps + 1 row for the initial state
		if( (newThin < x$thin) | (newThin %% x$thin) )
			stop(paste("thin.twDEMCPops: increased thin must be a positive multiple of former thin",x$thin))
		start <- max(start)	#
		if( start < 0)
			stop(paste("thin.twDEMCPops: argument start must be at least 0 but was ",start))
		nSPops <- getNSamples(x)
		thinFac <- newThin %/% x$thin
		startT <- ceiling( start / x$thin ) * x$thin		# adjust start time so that it coincides with next start of next thinning interval
		if( 1 == length(end) ) end <- rep(end, length(x$pops) )
		#nGenPops <- getNGen(x)
		#iPop=1
		for( iPop in seq_along(x$pops)){
			maxSampleTime <- iSample2time(nSPops[iPop], thin=x$thin)
			endi <- end[[iPop]]
			if( is.null(endi) || !is.finite(endi) || endi>maxSampleTime) endi=maxSampleTime 	
			if( endi < 1)
				stop(paste("thin.twDEMCPops: argument end must be at least 1 (one generation from 0 to 1) but was",endi))
			#thin own past: keep first line and every occurenc of multiple thin
			endT <- startT + floor( (endi-startT) / newThin) * newThin 			# adjust end time so that it coincides with beginning of thinning interval of end
			iStartEnd <- time2iSample( c(startT,endT), thin=x$thin, match="none" )
			iKeep <- seq(iStartEnd[1],iStartEnd[2],by=thinFac)
			#mtrace(subset.twDEMCPops)
			x <- tmp <- subset.twDEMCPops( x, iKeep, iPops=iPop )
		}
		x$thin <- newThin
		#time2iSample(70,5)
		#if( !is.null(x$nGenBurnin) ) res$nGenBurnin <- pmax(0,x$nGenBurnin-startT) 
		#if(!doKeepBatchCall) attr(res,"batchCall") <- NULL
		x
	})
#mtrace(thin.twDEMCPops)
attr(thin.twDEMCPops,"ex") <- function(){
	if( FALSE ){
		#mtrace(thin.twDEMCPops)
		tmp <- thin(res, start=4, newThin=8 )
		all.equals( c(1, 11), getNSamples(tmp))
	}
	data(twdemcEx1)
	x <- twdemcEx1
	c( nGen=getNGen(twdemcEx1), thin=twdemcEx1$thin, nSample=getNSamples(twdemcEx1), nGenBurnin=twdemcEx1$nGenBurnin )
	
	thinned <- thin(twdemcEx1, start=twdemcEx1$nGenBurnin)	# removing burnin period
	c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin )	#15 sample describing 70 generations
	
	thinned <- thin(twdemcEx1, start=twdemcEx1$nGenBurnin, newThin=10)	
	c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin )	#8 samples describing 70 generations
}


#setMethodS3("as.mcmc.list","twDEMCPops", function( 
as.mcmc.list.twDEMCPops <- function( 
	### Converts list of type twDEMCPops (result of \code{\link{twDEMCBlockInt}}) to coda's \code{mcmc.list}. 
	x				##<< the output of \code{\link{twDEMCBlockInt}}) run
	,...
	,useThinning=TRUE	##<< if TRUE thinning is used to make populations the same length, if FALSE they are cut to shortest population
	,minPopLength=NULL	##<< integer: if given, shorter populations are dropped 
){
	xStack <- concatPops.twDEMCPops(x, useThinning=useThinning, minPopLength=minPopLength)
	as.mcmc.list.twDEMC(xStack)
}


.unstackPopsTwDEMCPops <- function(
	### unstacks a stacked population (by \code{combineTwDEMCPops(mergeMethod="stack")})
	pop	##<< stacked population
	,nChain
){
	#print(".unstackPopsTwDEMCPops: start"); recover()
	# if nCases in pop is not a multiple of nChain, then last rows are skipped.
	nCases <- nrow(pop$parms) %/% nChain
	pops <- vector("list", nChain)
	#iChain <- nChain
	for( iChain in (1:nChain) ){
		#cases <- if( nCases==0) FALSE else nCases*(iChain-1) + (1:nCases)
		cases <- if( nCases==0) FALSE else (0:(nCases-1))*nChain+iChain	# all chains spread across initial chain
		pops[[iChain]] <- pop	# keeping entries spaceInd, upperParBound, lowerParBounds, and splits		
		pops[[iChain]]$parms <- pop$parms[cases,, ,drop=FALSE] 
		pops[[iChain]]$logDen <- pop$logDen[cases,, ,drop=FALSE] 
		pops[[iChain]]$resLogDen <- pop$resLogDen[cases,, ,drop=FALSE] 
		pops[[iChain]]$pAccept <- pop$pAccept[cases,, ,drop=FALSE] 
		pops[[iChain]]$temp <- pop$temp[cases, ,drop=FALSE] 
	}
	pops
	### list of subsets (by rows)
}

.concatChainsTwDEMCPops <- function(
	### combine given populations to one big population with more chains  
	pops
){
	newPop <- pops[[1]]	# keeping entries spaceInd
	if( length(pops) > 1){
		along = 3 
		newPop$parms <- abind( lapply(pops,"[[","parms"), along=along ) 
		newPop$logDen <- abind( lapply(pops,"[[","logDen"), along=along )
		newPop$resLogDen <- abind( lapply(pops,"[[","resLogDen"), along=along )
		newPop$pAccept <- abind( lapply(pops,"[[","pAccept"), along=along )
		#only one temp for all chain and pops newPop$temp <- abind( lapply(pops,"[[","temp"), along=2 ) 
		
		#---- update the parameter bounds
		pB <- .parBoundsEnvelope(pops)
		newPop[ names(pB)] <- pB	
		##details<<
		## entry splits is discarded, as it is not generally determined.
		newPop$splits <- numeric(0)	 
		newPop
	}
	newPop
}

.mergePopTwDEMC <- function(
	### merge population iPop to other populations
	pops	##<< list of populations with entries upperParBounds, lowerParBounds, spaceInd, splits
	,iPop	##<< index of the population that should be merged to other ones
	,pPops	##<< the initial proabilities of the subspaces
	,mergeMethod="random"	##<< the mergeMethod of \code{\link{combineTwDEMCPops}}
){
	#print(".mergePopTwDEMC"); recover()
	#stop(".mergePopTwDEMC: not implemented yet.")
	if( length(pops) <= 1 ){
		warning(".mergePopTwDEMC: called with only one population.")
		return( list(pops=pops, pPops=pPops) )
	}
	iPop0 <- iPop	# index of population in all populations of all space replicates
	popSource <- pops[[iPop]]	# the population to merge from
	spacesPop <- tmp <- sapply( pops, "[[", "spaceInd" )
	iSameSpace <- which(spacesPop == popSource$spaceInd)	
	iOtherSpace <- which(spacesPop != popSource$spaceInd)
	newPops <- pops[ iSameSpace]	# possible populations to merge to of same space replicate
	newPPops <- pPops[iSameSpace]
	# translate iPop0 in all populations into iPop index within same space 
	nPop <- length(newPops)
	tmp[iSameSpace] <- 1:nPop
	tmp[-iSameSpace] <- NA	# for detecting errors
	iSPop <- tmp[iPop0]		# index of iPop in newPops (of same space)
	#determine neighbouring pops with same first part of splits
	#pop <- popsS[-iPop][1]
	jPops <- (1:nPop)[-iSPop][ sapply( newPops[-iSPop], function(pop){ 
			(length(pop$splits) >= length(popSource$splits)) && all(pop$splits[ 1:length(popSource$splits)] == popSource$splits)
		})]
	# j = jPops[1]
	if( 1 != length(jPops) ){
		# need to split the population to merge
		# in addition, some of the spaces might have no common border (seperated by a further split)
		# update newPops[jPops] and also strip non-border indices from jPops
		jPopsBorder <- integer(max(jPops))		# initialized by 0, not updated means no border
		#lapply( newPops, "[[", "splits" )
		#.getParBoundsPops(c(subSpacesAll,list(popM)))
		# if upper Bound of source pop equals lower bound of target pop, decrease lower target bound 
		#iB=1
		for( iB in seq_along(popSource$upperParBounds) ){
			parName <- names(popSource$upperParBounds)[iB]
			ubVal <- popSource$upperParBounds[iB]
			lbVal <- if( 0 != length(popSource$lowerParBounds) && is.finite(popSource$lowerParBounds[parName])) 
				popSource$lowerParBounds[parName] else NULL
			#jPop = jPops[1]
			#sapply(newPops[jPops], "[[", "splits")
			for( jPop in jPops ){
				lb <- newPops[[jPop]]$lowerParBounds
				if( 0 != length(lb) && is.finite(lb[parName]) && lb[parName]==ubVal ){
					if( is.null(lbVal) )
						newPops[[jPop]]$lowerParBounds <- newPops[[jPop]]$lowerParBounds[-match(parName,names(newPops[[jPop]]$lowerParBounds))]					
					else
						newPops[[jPop]]$lowerParBounds[parName] <- lbVal					
					jPopsBorder[jPop] <- jPop
				}
			}
		}
		# if lower Bound of source pop equals upper bound of target pop, increase upper target bound 
		for( iB in seq_along(popSource$lowerParBounds) ){
			parName <- names(popSource$lowerParBounds)[iB]
			lbVal <- popSource$lowerParBounds[iB]
			ubVal <- if( 0 != length(popSource$upperParBounds) && is.finite(popSource$upperParBounds[parName])) 
					popSource$upperParBounds[parName] else NULL
			#jPop = 1
			for( jPop in jPops ){
				ub <- newPops[[jPop]]$upperParBounds
				if( 0 != length(ub) && is.finite(ub[parName]) && ub[parName]==lbVal ){
					if( is.null(ubVal) )
						newPops[[jPop]]$upperParBounds <- newPops[[jPop]]$upperParBounds[-match(parName,names(newPops[[jPop]]$upperParBounds))]					
					else
						newPops[[jPop]]$upperParBounds[parName] <- ubVal					
					jPopsBorder[jPop] <- jPop
				}
			}
		}
		#mtrace(divideTwDEMCPop)
		jPopsSplitBorder <- jPopsBorder[ jPopsBorder != 0]
		jPopsSplitNonBorder <- setdiff( jPops, jPopsSplitBorder )	# remember those non neighbouring to remove the split
	} else{
		# nPops == 1
		jPopsSplitBorder <- jPops
		jPopsSplitNonBorder <- integer(0) 
	} 
	#	
	# jPops indicates position of subspaces in newPops, i.e those populations that share a common splitting history with the merge pop
	# jPopsSplitBorder is the subset of jPops that have a common border
	# jPopsSplitNonBorder is complementing subset of jPops that have no common border
	#
	pSubsSource <- 1
	subsPopSource <- if( length(jPopsSplitBorder) > 1 ){
		# must merge this population to several other populations, need to split before
		subsPopSource <- divideTwDEMCPop( popSource, newPops[jPopsSplitBorder], isCheckBorderConsistency=FALSE ) #popSource has more constrained support, hence do not check borders
		# may loose a few samples on dividing into subspaces because all chains must have equal length
		#if( nrow(popSource$parms) != sum( sapply( subsPopSource,function(pop){ nrow(pop$parms)})) )
		#	stop(".mergePopTwDEMC: sum of samples of subspaces does not match the original sample number.")
		.nS <- sapply(subsPopSource, function(pop){ nrow(pop$parms) })
		pSubsSource <- .nS/sum(.nS)
		subsPopSource
	}else{
		list(popSource)
	}
	if( all(is.na(pSubsSource)) ) # splitting a very small pop, may loose all samples, rendergin pSubsSource NA (division by zero)
		pSubsSource[] <- 1/length(pSubsSource)		
	#nrow(popSource$parms)
	#sapply( lapply(subsPopSource,"[[","parms" ), nrow )
	#
	#-------- do the actual merge
	#ij = 1
	for( ij in seq_along(jPopsSplitBorder) ){
		jPop <- jPopsSplitBorder[ij]
		popDest <- newPops[[jPop]]
		subPopSource <- subsPopSource[[ij]]
		# c(popDest$lowerParBounds, popDest$upperParBounds )
		# c(subPopSource$lowerParBounds, subPopSource$upperParBounds )  # a combination of uB and lB should match
		# mtrace(combineTwDEMCPops)
		if( nrow(subPopSource$parms) != 0){
			resCombine <-  combineTwDEMCPops( list(popDest, subPopSource), mergeMethod=mergeMethod )
			newPops[[jPop]] <-	if( resCombine$popCases[ (.nC <- length(resCombine$popCases))] == 2){
				# make sure that last entry is from popDest (first pop in combineTwDEMCPops) to continue with current state
				# so switch the last entry from popDest with last entry of the combined pop
				iSwitch <- .nC+1- match( 1, rev(resCombine$popCases) ) # position of the last 1
				.tmp1 <- .subsetTwDEMCPop( resCombine$pop, iKeep=seq_along(resCombine$popCases)[-iSwitch] )
				.tmp2 <- .subsetTwDEMCPop( resCombine$pop, iKeep=seq_along(resCombine$popCases)[iSwitch] )
				.tmp <- combineTwDEMCPops( list(.tmp1, .tmp2), mergeMethod="stack" )$pop  			
			}else .tmp <- resCombine$pop
	if( nrow(newPops[[jPop]]$parms) != (nrow(popDest$parms) + nrow(subPopSource$parms))) stop(".mergePopTwDEMC: encountered mismatch in sample numbers.")	
			#.getParBoundsPops(c( list(subPopSource), list(popDest), list(newPops[[jPop]])) )
			# c(tmp$lowerParBounds, tmp$upperParBounds )  # matching border should disappear
		}
		newPPops[jPop] <- newPPops[jPop] + newPPops[iSPop]*pSubsSource[ij]   # update the probability of the space
		#c( nrow(pop$parms), nrow(popM$parms), nrow(tmp$parms) )
		newPops[[jPop]]$splits <- popDest$splits[ -length(popSource$splits) ]  # we just removed this splitting point
	}
	#
	#-------- remove the split from all non-bordering populations with same split history too
	#rmSplit <-  popSource$splits[ length(popSource$splits) ]
	#jPop <- 7
	#jPop <- jPopsSplitNonBorder[1]
	# sapply( newPops[jPopsSplitNonBorder], "[[", "splits" )
	for( jPop in jPopsSplitNonBorder ){
		newPops[[jPop]]$splits <- newPops[[jPop]]$splits[ -length(popSource$splits) ]  
	}
	#
	# ------- delete the merged source population 
	# only after merging so that indices hold true
	newPops[[iSPop]] <- NULL	
	newPPops <- newPPops[-iSPop]
	iSameSpace0 <- iSameSpace
	iSameSpace <- iSameSpace[ iSameSpace != iPop0 ]
	# newPops has changed, update the position indices
	bo <- jPops > iSPop; jPops[bo] <- jPops[bo]-1 
	bo <- jPopsSplitBorder > iSPop; jPopsSplitBorder[bo] <- jPopsSplitBorder[bo]-1 
	bo <- jPopsSplitNonBorder > iSPop; jPopsSplitNonBorder[bo] <- jPopsSplitNonBorder[bo]-1 
	#
	#-------- check consistency
	# may have lost some samples during splitting because all chains in one pop need to have the same lenght
	.nS0 <- sapply( pops[ iSameSpace0 ], function(pop){ nrow(pop$parms)})
	.nS <- sapply( newPops, function(pop){ nrow(pop$parms)})
	if( sum(.nS) > sum(.nS0)  ){
		stop(".mergePopTwDEMC: sample number in merged populations is larger than original.")
		nrow( popSource$parms )
		sapply( lapply(subsPopSource,"[[","parms" ), nrow )
	}
	if( !isTRUE(all.equal( sum(newPPops),1)) )
		stop(".mergePopTwDEMC: probabilities of space with merged does not sum to 1")
	#
	##value<< list with entries
	resMerge <- list(
		##describe<<
		pops = c(pops[iOtherSpace], newPops)        ##<< list (nPops): the pops with one population less, that has been merged to the other pops
		,pPops = c( pPops[iOtherSpace], newPPops )	##<< numeric vector (nPops): the updated proportions of the populations within space
		#,nSub = length(subsPopSource)				##<< the number of subpopulation thats the population was splitted into
		,iPopsBefore = c(iOtherSpace, iSameSpace) 	##<< integer vector (nPops): index of the population in original set of populations
		,iPopsModified = length(iOtherSpace) + seq_along(newPops)[jPopsSplitBorder]	##<< integer vector (nPopSameBorder): index of the populations that comprise a part of the merged population
		,pPopsSource = pSubsSource					##<< numeric vector (nPopSameBorder): proportion of the merged population that is part of the new population. Indices correspond to iPopsNew 
		)	
		##end<<
	# lapply(pops[-iSameSpace], "[[", "splits")
	# sapply(pops[-iSameSpace], "[[", "spaceInd")
	# lapply(newPops, "[[", "splits")
	# sapply(newPops, "[[", "spaceInd")
}
attr(.mergePopTwDEMC,"ex") <- function(){
	data(den2dCorEx)
	#mtrace(divideTwDEMCSteps)
	#trace("twDEMCBlockInt", recover)
	mc0 <- den2dCorEx$mcSubspaces0
	(tmp <- lapply(mc0$pops, "[[", "splits") )
	pops <- mc0$pops
	#spaceI <- den2dCorEx$subspaces0[[1]] 
	pPops <- do.call( c, lapply( den2dCorEx$subspaces0, function(spaceI){
			sapply( spaceI$spaces, "[[", "pSub")
		}))
	getSpacesPop(mc0)
	.getParBoundsPops(mc0$pops)

	iPop <- 8	# second space replicate, must merge to 9 only
	#mtrace(.mergePopTwDEMC)
	resMerge <- .mergePopTwDEMC( mc0$pops, iPop=iPop, pPops=pPops)
	all.equal(13, length(resMerge$pPops) )
	all.equal(13, length(resMerge$pops) )
	all.equal(nrow(pops[[8]]$parms)+nrow(pops[[9]]$parms), nrow(resMerge$pops[[8]]$parms) ) 
	all.equal(pPops[[8]]+pPops[[9]], resMerge$pPops[[8]] ) 

	iPop <- 10	# second space replicate, hast split with 8 and 9 but must merge only to second
	#mtrace(.mergePopTwDEMC)
	resMerge <- .mergePopTwDEMC( mc0$pops, iPop=iPop, pPops=pPops)
	all.equal(13, length(resMerge$pPops) )
	all.equal(13, length(resMerge$pops) )
	all.equal(nrow(pops[[10]]$parms)+nrow(pops[[9]]$parms), nrow(resMerge$pops[[9]]$parms) ) 
	all.equal(pPops[[10]]+pPops[[9]], resMerge$pPops[[9]] )
	#sapply( pops, function(pop){ nrow(pop$parms) })
	#sapply( resMerge$pops, function(pop){ nrow(pop$parms) })
	
	iPop <- 3	# from first space replicate, must merge to 1 and 2 which share the same upper bound on a
	#mtrace(.mergePopTwDEMC)
	resMerge <- .mergePopTwDEMC( mc0$pops, iPop=iPop, pPops=pPops)
	all.equal(13, length(resMerge$pPops) )
	all.equal(13, length(resMerge$pops) )
	all.equal(nrow(pops[[1]]$parms)+nrow(pops[[2]]$parms)+nrow(pops[[3]]$parms), nrow(resMerge$pops[[8]]$parms)+nrow(resMerge$pops[[9]]$parms) ) 
	all.equal(pPops[[1]]+pPops[[2]]+pPops[[3]], resMerge$pPops[[8]]+resMerge$pPops[[9]] )
	pops[[1]]$splits
	.getParBoundsPops(resMerge$pops[8:9])
}


combineTwDEMCPops <- function(
	### combine given populations to one big chain with more cases
	pops					##<< list of population objects
	,popCases=integer(0)	##<< integer vector (sum(getNSamples(pops))): specifying for each case (row) from which population it is filled
	,mergeMethod="random"	##<< method of merging the pops, see \code{\link{twMergeSequences}}
	,nInSlice=4		##<< sequence length from each population (only for mergeMethod \code{slice}
){
	#print("start of combineTwDEMCPops"); recover()
	nPop <- length(pops)
	newPop <- pop1 <- pops[[1]]
	if( nPop > 1){
		nSamplePop <- sapply(pops,function(pop){nrow(pop$parms)})
		nSample <- sum(nSamplePop)
		if( 0 == length(popCases))
			popCases <- twMergeSequences(nSamplePop, mergeMethod=mergeMethod, nInSlice=nInSlice )
		if( nSample != length(popCases) )
			stop("combineTwDEMCPops: length of argument popCases must equal sum of samples of all populations")
		#----  initialize entries
		tmp <- "parms"; newPop[[tmp]] <- array( NA_real_, dim=c(nSample,dim(pop1[[tmp]])[2:3]), dimnames=dimnames(pop1[[tmp]])) 
		tmp <- "logDen"; newPop[[tmp]] <- array( NA_real_, dim=c(nSample,dim(pop1[[tmp]])[2:3]), dimnames=dimnames(pop1[[tmp]])) 
		tmp <- "resLogDen"; newPop[[tmp]] <- array( NA_real_, dim=c(nSample,dim(pop1[[tmp]])[2:3]), dimnames=dimnames(pop1[[tmp]])) 
		tmp <- "pAccept"; newPop[[tmp]] <- array( NA_real_, dim=c(nSample,dim(pop1[[tmp]])[2:3]), dimnames=dimnames(pop1[[tmp]])) 
		tmp <- "temp"; newPop[[tmp]] <- matrix( NA_real_, nrow=nSample, ncol=ncol(pop1[[tmp]]),dimnames=dimnames(pop1[[tmp]])) 
		#---- copy the entries
		#iPop=1
		for( iPop in 1:nPop){
			pop <- pops[[iPop]]
			is <- which(popCases==iPop)
			tmp <- "parms"; newPop[[tmp]][is,,] <- pop[[tmp]]
			tmp <- "logDen"; newPop[[tmp]][is,,] <- pop[[tmp]]  
			tmp <- "resLogDen"; newPop[[tmp]][is,,] <- pop[[tmp]]  
			tmp <- "pAccept"; newPop[[tmp]][is,,] <- pop[[tmp]]  
			tmp <- "temp"; newPop[[tmp]][is,] <- pop[[tmp]]  
		} # iPop
		# newPop$parms[,1,1]
		
		#---- update the parameter bounds
		# .getParBoundsPops(pops)
		# mtrace(.parBoundsEnvelope)
		pB <- .parBoundsEnvelope(pops)
		newPop[ names(pB)] <- pB	
		##details<<
		## entry splits is discarded, as it is not generally determined.
		newPop$splits <- numeric(0)	 
	} # length(pops) > 1
	##value<<
	list( 
		pop=newPop			##<< the merged population
		,popCases=popCases 	##<< integer vector (nSample): population that case is taken from
	)
}
attr(combineTwDEMCPops,"ex") <- function(){
	data(den2dCorEx)
	runs0 <- subPops(den2dCorEx$mcSubspaces0, iSpaces=c(1))
	getSpacesPop(runs0)	# all populatins belonging to space replicate 1
	#pops <- runs0$pops
	#mtrace(combineTwDEMCPops)
	res <- combineTwDEMCPops(runs0$pops) 
	#str(res)
	head(res$pop$parms[,,1])
	
	res <- combineTwDEMCPops(runs0$pops[3:4])
	#str(res)
	
}

setMethodS3("stackPopsInSpace","twDEMCPops", function( 
		### Combine populations for subspaces to bigger populations 
		x
		,...	##<< arguments passed \code{\link{combineTwDEMCPops}} such as \code{mergeMethod=stack/slice/random}
		,spacePop = getSpacesPop(x)
	){
		# stackPopsInSpace.twDEMCPops
		##seealso<<   
		## \code{\link{stackPops.twDEMCPops}}
		## \code{\link{combineTwDEMCPops}}
		## \code{\link{stackChains.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		#iSpace=1
		resComb <- lapply( unique(spacePop), function(iSpace){
			iPops <- which(spacePop == iSpace)
			combineTwDEMCPops(x$pops[iPops], ...) 
		})
		x$pops <- lapply( resComb, "[[", "pop" )
		### twDEMCPops object with nSpace populations
		x
	})
attr(stackPopsInSpace.twDEMCPops,"ex") <- function(){
	data(den2dCorEx)
	getNSamples(den2dCorEx$mcBulk)
	res <- stackPopsInSpace( den2dCorEx$mcSubspaces0 )
	getNSamples(res)	# lost a few samples in sorting chains to subspaces
	#mtrace(concatPops.twDEMCPops)
	plot( as.mcmc.list(den2dCorEx$mcBulk), smooth=FALSE ) # original before splitting into subspaces
	plot( as.mcmc.list(res), smooth=FALSE )		# stacked populations of subspaces
	plot( as.mcmc.list(stackPopsInSpace( den2dCorEx$mcSubspaces0,mergeMethod="stack" )), smooth=FALSE )
	plot( as.mcmc.list(stackPopsInSpace( den2dCorEx$mcSubspaces0,mergeMethod="slice" )), smooth=FALSE )
}

setMethodS3("stackPops","twDEMCPops", function( 
		### Combine all populations (across subspaces) to one big population 
		x
		,...	##<< arguments passed \code{\link{combineTwDEMCPops}} such as \code{mergeMethod=stack/slice/random}
	){
		# stackPopsInSpace.twDEMCPops
		##seealso<<   
		## \code{\link{stackPopsInSpace.twDEMCPops}}
		## \code{\link{combineTwDEMCPops}}
		## \code{\link{stackChains.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		#iSpace=1
		x$pops <- list(combineTwDEMCPops(x$pops, ...)$pop)
		x
	})


#------------------------ 
.stackChainsWithinPop <- function(
	pop
	,mergeMethod="stack"
){
	nc <- dim(pop$parms)[3]
	nr <- nrow(pop$parms)
	newPop <- pop
	if( nc != 1 ){
		newPop$parms <- abind( lapply( 1:nc, function(iChain){ pop$parms[,,iChain ,drop=FALSE]}), along=1 ) 
		newPop$logDen <- abind( lapply( 1:nc, function(iChain){ pop$logDen[,,iChain ,drop=FALSE]}), along=1 )
		newPop$resLogDen <-abind( lapply( 1:nc, function(iChain){ pop$resLogDen[,,iChain ,drop=FALSE]}), along=1 )
		newPop$pAccept <- abind( lapply( 1:nc, function(iChain){ pop$pAccept[,,iChain ,drop=FALSE]}), along=1 )
		newPop$temp <- abind( lapply( 1:nc, function(iChain){ pop$temp }), along=1 )  # repeat the temperature
		newPop$Y <- abind( lapply( 1:nc, function(iChain){ pop$Y[,,iChain ,drop=FALSE]}), along=1 )
		if( mergeMethod != "stack" && nr != 0){
			# need to resort the entries
			chainInd <- twMergeSequences( rep(nr,nc), mergeMethod=mergeMethod )
			chainIndY <- twMergeSequences( rep(nrow(pop$Y),nc), mergeMethod=mergeMethod )
			#iChain <- 1
			for( iChain in 1:nc){
				iPos <- which(chainInd==iChain)
				iPosY <- which(chainIndY==iChain)
				newPop$parms[iPos,,1] <- pop$parms[,,iChain]
				newPop$logDen[iPos,,1] <- pop$logDen[,,iChain]
				newPop$resLogDen[iPos,,1] <- pop$resLogDen[,,iChain]
				newPop$pAccept[iPos,,1] <- pop$pAccept[,,iChain]
				newPop$temp[iPos,] <- pop$temp
				newPop$Y[iPosY,,1] <- pop$Y[,,iChain]
			}
		}
	}
	newPop
}

setMethodS3("stackChainsPop","twDEMCPops", function( 
		### Combine MarkovChains of each population of a twDEMC to a single chain.
		x
		,mergeMethod="stack"
		,...
	){
		#stackChainsPop.twDEMCPops
		##seealso<<   
		# stack logDen for each population
		#nPop = getNPops(x)
		#iPop = nPop
		#pop <- x$pops[[nPop]]
		xNew <- x
		xNew$pops <- newPops <- lapply( x$pops, .stackChainsWithinPop, mergeMethod=mergeMethod )
		xNew
	})
attr(stackChainsPop.twDEMCPops,"ex") <- function(){
	data(twdemcEx1)
	getNChainsPop(twdemcEx1)	# four chains within population
	getNSamples(twdemcEx1)		# with 26 samples
	
	res <- stackChainsPop(twdemcEx1)
	getNChainsPop(res)			# only 1 chains
	getNSamples(res)			# but with 26*4=104 samples
	str(concatPops(res))
}


sumLogDenComp <- function(
	### sum LogDen components within density
	resLogDen		##<< numeric matrix (nStep x nResComp x nChain): of log-Densities
	,dInfos=NULL	##<< list of densities each a list with entry resCompPos: integer vector, specifying the result components within density
		##<< If dInfos is of length 0, all components are summed and a matrix instead of a array is returned
){
	if( 0 == length(dInfos) ){
		# sum sum within chain
		sumD <- apply(resLogDen, c(1,3), sum )
		dimnames(sumD) <- dimnames(resLogDen)[c(1,3)]
		sumD
	}else{
		#iInfo <- 1
		logDenSumL <- lapply( seq_along(dInfos), function(iInfo){
			resCompPos <- dInfos[[iInfo]]$resCompPos
			if( 0 == length(resCompPos)) stop("sumLogDenComp: each entry of dInfos must provide entry resCompPos of length > 0")
			if( 1 == length(resCompPos)) adrop( resLogDen[,resCompPos, ,drop=FALSE],2) 
			tmp <- apply(resLogDen[,resCompPos, ,drop=FALSE], c(1,3), sum )
		})
		logDenSumArr1 <- abind(logDenSumL, rev.along=0)
		dimnames(logDenSumArr1)[3] <- names(dInfos)
		sumD <- aperm(logDenSumArr1, c(1,3,2))
	}
	##value<< numeric array (nStep x nDensities x nChain )
	##seealso<< \code{\link{calcTemperatedLogDenChains.array}}
}





