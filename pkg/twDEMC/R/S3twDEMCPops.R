

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
		## \code{\link{getNChains.twDEMCPops}}
		## \code{\link{getNChainsPop.twDEMCPops}}
		## \code{\link{getNParms.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCInt}}
		#mtrace(getNSamples.twDEMCPops)
		(getNSamples(x)-1)*x$thin
		### integer vector, number of completed generations
	})
attr( getNGen.twDEMCPops, "ex") <- function(){
	data(twdemcEx1)
	getNGen(twdemcEx1)
	getNSamples(twdemcEx1)
	twdemcEx1$thin
	getNPops(twdemcEx1)
	getNChains(twdemcEx1)
	getNChainsPop(twdemcEx1)
	getNParms(twdemcEx1)
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
		## ,\code{\link{twDEMCInt}}
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
		## ,\code{\link{twDEMCInt}}
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
		## ,\code{\link{twDEMCInt}}
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
		## ,\code{\link{twDEMCInt}}
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
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subset.twDEMCPops}}
		## ,\code{\link{twDEMCInt}}
		##details<< There is only one sample per thinning interval of length \code{x$thin}.
		sapply( x$pops, function(pop){ nrow(pop$parms) })
		### integer vector: number of samples in each population of twDEMCPops
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
	## ,\code{\link{twDEMCInt}}
	nStepsPop <- getNSamples(x)
	if( 1 == length(minPopLength) ){
		iKeep <- which( nStepsPop >= minPopLength )
		x <- subPops(x, iPops=iKeep )
		nStepsPop <- nStepsPop[iKeep]
	}
	nSteps <- min(nStepsPop)
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
	x$temp <- structure( abind( lapply(pops,"[[","temp"), along=2), dimnames=list(steps=NULL,pops=NULL) )
	x$pAccept <- structure( abind( lapply(pops,"[[","pAccept"), along=3), dimnames=dimnames(p1$pAccept))
	x$resLogDen <- structure( abind( lapply(pops,"[[","resLogDen"), along=3), dimnames=dimnames(p1$resLogDen))
	x$logDen <- structure( abind( lapply(pops,"[[","logDen"), along=3), dimnames=dimnames(p1$logDen))
	YL <-  lapply(pops,"[[","Y")
	nY <- min(sapply(YL,nrow))
	YLs <- lapply(YL, function(Y){ Y[nrow(Y)-((nY-1):0),,,drop=FALSE] })
	x$Y <- structure( abind(YLs, along=3), dimnames=dimnames(p1$Y))
	x$upperParBoundsPop = lapply( pops, "[[", "upperParBounds" )
	x$lowerParBoundsPop = lapply( pops, "[[", "lowerParBounds" )
	class(x) <- c("list","twDEMC")
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

#mtrace(concatPops.twDEMCPops)

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
	## \code{\link{twDEMCInt}}
	
	##details<< 
	## There are several methods access properties a result of an \code{\link{twDEMCInt}}, i.e.
	## an object of class \code{twDEMCPops} \itemize{
	## \item{ number of generations: \code{\link{getNGen.twDEMCPops}}  } 
	## \item{ number of samples (only one sample each thinning inteval): \code{\link{getNSamples.twDEMCPops}}  } 
	## \item{ number of chains: \code{\link{getNChains.twDEMCPops}}  } 
	## \item{ number of populations: \code{\link{getNPops.twDEMCPops}}  } 
	## \item{ number of chains per population: \code{\link{getNChainsPop.twDEMCPops}}  } 
	## \item{ number of parameters: \code{\link{getNParms.twDEMCPops}}  } 
	## \item{ thinning interval: \code{res$thin}  } 
	##}
	##
	## There are several methods to transform or subset the results of an \code{\link{twDEMCInt}} run. \itemize{
	## \item{ select chains or sub-populations: this method  } 
	## \item{ thin all the chains: \code{\link{thin.twDEMCPops}}  } 
	## \item{ select subset of cases: \code{\link{subset.twDEMCPops}}  }
	## \item{ combine several twDEMCPops results to a bigger set of populations \code{\link{combinePops.twDEMCPops}}  }
	## \item{ stack all the results of all chains to one big matrix \code{\link{stackChains.twDEMCPops}}  } 
	##}
	##
	## There are several methods utilize the functions of the coda package. \itemize{
	## \item{ convert an twDEMCPops to a coda mcmc.list \code{\link{as.mcmc.list.twDEMCPops}}  } 
	## \item{ applying a function to all of the chains: \code{\link{mcmcListApply}}  }
	## \item{ stack all the results of all chains to one big matrix: \code{\link{stackChains.mcmc.list}}  } 
	## \item{ plotting a subset of the chains and cases: \code{\link{plotThinned.mcmc.list}}  } 
	## \item{ transforming parameters \code{\link{transOrigPopt.mcmc.list}}  } 
	##}
	
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
		x$pops[[iPop]]$parms <- x$pops[[iPop]]$parms[iKeep,,, drop=FALSE] 
		x$pops[[iPop]]$logDen <- x$pops[[iPop]]$logDen[iKeep,,, drop=FALSE] 
		x$pops[[iPop]]$resLogDen <- x$pops[[iPop]]$resLogDen[iKeep,,, drop=FALSE] 
		x$pops[[iPop]]$pAccept <- x$pops[[iPop]]$pAccept[iKeep,,, drop=FALSE]
		x$pops[[iPop]]$temp <- x$pops[[iPop]]$temp[iKeep, drop=FALSE]
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

setMethodS3("subPops","twDEMCPops", function( 
		### Condenses an twDEMCPops List to the chains iChains e.g \code{1:4}.
		x
		, iPops 		##<< populations to keep	 
		,... 
	#, doKeepBatchCall=FALSE	##<< wheter to retain the batch call attribute of x
	){
		##seealso<<   
		## \code{\link{twDEMCInt}}
		## \code{\link{subset.twDEMCPops}}
		x$pops <- x$pops[iPops]
		x
	})
#mtrace(subPops.twDEMCPops)



setMethodS3("squeeze","twDEMCPops", function(
	### Reduces the rows of populations so that all populations have the same number of samples. 
	x, ##<< the twDEMCPops list to thin 
	...,
	length.out=min(getNSamples(x)),	##<< number of steps in each population
	dropShortPops=FALSE				##<< if set to TRUE, pops with less samples than length.out are dropped
){
	# squeeze.twDEMCPops
	nSamplesPop <- getNSamples(x)
	if( dropShortPops ){
		x <- subPops( x, iPops=which(nSamplesPop < length.out))
	}else if( length.out > min(nSamplesPop)) stop(
			"squeeze.twDEMCPops: specified a length that is longer than the shortest population. Use dropShortPops=TRUE to drop shorter populations.")
	##details<< 
	## all populations with 
	for( iPop in seq_along(x$pops) ){
		iKeep <- seq(1,nSamplesPop[iPop],length.out=length.out) 
		x$pops[[iPop]]$parms <- x$pops[[iPop]]$parms[iKeep,,, drop=FALSE] 
		x$pops[[iPop]]$logDen <- x$pops[[iPop]]$logDen[iKeep,,, drop=FALSE] 
		x$pops[[iPop]]$resLogDen <- x$pops[[iPop]]$resLogDen[iKeep,,, drop=FALSE] 
		x$pops[[iPop]]$pAccept <- x$pops[[iPop]]$pAccept[iKeep,,, drop=FALSE]
		x$pops[[iPop]]$temp <- x$pops[[iPop]]$temp[iKeep, drop=FALSE]
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

setMethodS3("thin","twDEMCPops", function( 
		### Reduces the rows of an twDEMCPops object (list returned by \code{\link{twDEMCInt}}) to correspond to a thinning of \code{newThin}.
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
			x <- subset.twDEMCPops( x, iKeep, iPops=iPop )
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
	### Converts list of type twDEMCPops (result of \code{\link{twDEMCPops}}) to coda's \code{mcmc.list}. 
	x				##<< the output of \code{\link{twDEMCBlockInt}}) run
	,...
	,useThinning=TRUE	##<< if TRUE thinning is used to make populations the same length, if FALSE they are cut to shortest population
	,minPopLength=NULL	##<< integer: if given, shorter populations are dropped 
){
	xStack <- concatPops.twDEMCPops(x, useThinning=useThinning, minPopLength=minPopLength)
	as.mcmc.list.twDEMC(xStack)
}

