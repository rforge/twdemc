#	(tres <- twUtestF("S3twDEMC")) #thin, subset, subChains, combinePops

R.methodsS3::setMethodS3("subChains","twDEMC", function( 
		### Condenses an twDEMC List to the chains iChains e.g \code{1:4}.
		x, iChains=NULL, iPops=NULL, 
		nPop=NULL, ##<< number of populations in x
		... 
		, doKeepBatchCall=FALSE	##<< wheter to retain the batch call attribute of x
){
		##seealso<<   
		## \code{\link{twDEMC}}
		
		##details<< 
		## There are several methods access properties a result of an \code{\link{twDEMCBlockInt}}, i.e.
		## an object of class \code{twDEMC} \itemize{
		## \item{ number of generations: \code{\link{getNGen.twDEMC}}  } 
		## \item{ number of samples (only one sample each thinning inteval): \code{\link{getNSamples.twDEMC}}  } 
		## \item{ number of chains: \code{\link{getNChains.twDEMC}}  } 
		## \item{ number of populations: \code{\link{getNPops.twDEMC}}  } 
		## \item{ number of chains per population: \code{\link{getNChainsPop.twDEMC}}  } 
		## \item{ number of parameters: \code{\link{getNParms.twDEMC}}  } 
		## \item{ thinning interval: \code{res$thin}  } 
		##}
		##
		## There are several methods to transform or subset the results of an \code{\link{twDEMCBlockInt}} run. \itemize{
		## \item{ select chains or sub-populations: this method  } 
		## \item{ thin all the chains: \code{\link{thin.twDEMC}}  } 
		## \item{ select subset of cases: \code{\link{subset.twDEMC}}  }
		## \item{ combine several twDEMC results to a bigger set of populations \code{\link{combinePops.twDEMC}}  }
		## \item{ stack all the results of all chains to one big matrix \code{\link{stackChains.twDEMC}}  } 
		##}
		##
		## There are several methods utilize the functions of the coda package. \itemize{
		## \item{ convert an twDEMC to a coda mcmc.list \code{\link{as.mcmc.list.twDEMC}}  } 
		## \item{ applying a function to all of the chains: \code{\link{mcmcListApply}}  }
		## \item{ stack all the results of all chains to one big matrix: \code{\link{stackChains.mcmc.list}}  } 
		## \item{ transforming parameters \code{\link{transOrigPopt.mcmc.list}}  } 
		##}
		#? 		## \item{ plotting a subset of the chains and cases: \code{\link{plotThinned.mcmc.list}}  } 


		# To help re-initializing the arguments to fLogDen \itemize{
		# \item{ transforming parameters \code{\link{subsetArgsFLogDen}}  }}
		#
		
		##details<< 
		## Alternatively to specification of iChains, one can specify a vector of populations and the total number of populations.
		if( is.null(nPop) ) nPop=getNPops(x)
		nChainPop = getNChainsPop(x)
		res <- x
		if( !is.null(iPops)){
			iChains = unlist( lapply( iPops, function(iPop){ (iPop-1)*nChainPop + (1:nChainPop) }))
			#same across all populations   res$temp <- x$temp[,iPops ,drop=FALSE]	#by nPop
			if( !is.null(x$nGenBurnin) ) res$nGenBurnin <- x$nGenBurnin[iPops]
		}else{
			# no populatios given: reduce to one population
			#res$temp <- matrix( x$temp[,1], nrow=nrow(x$temp), ncol=1, dimnames=dimnames(x$temp) )
			if( !is.null(x$nGenBurnin) ) res$nGenBurnin <- max(x$nGenBurnin)
		}
		res$parms <- x$parms[,,iChains, drop=FALSE]
		res$logDenComp <- x$logDenComp[,,iChains, drop=FALSE]
		res$logDen <- x$logDen[,,iChains, drop=FALSE] 
		res$pAccept <- x$pAccept[,,iChains, drop=FALSE]
		res$thin <- x$thin
		if( 0 < length(x$Y))
			res$Y <- x$Y[,,iChains, drop=FALSE]
		res$thin <- x$thin
		if(!doKeepBatchCall) attr(res,"batchCall") <- NULL
		res
		### a list of class twDEMC (see \code{\link{twDEMCBlockInt}})
	})
#mtrace(subChains.twDEMC)

iSample2time <- function(
	### Convert sample number to time given the thinning interval
	iSample	##<< integer vector: indices of the samples in recorded states
	,thin=1	##<< integer scalar: thinning interval
){
	##seealso<<   
	## \code{\link{time2iSample}}
	## \code{\link{getNGen.twDEMC}}
	## ,\code{\link{subChains.twDEMC}}
	## ,\code{\link{twDEMCBlockInt}}
	
	##details<<
	## Sample 1 corresponds to time zero
	## After first generation, the time sample 2 corresponds to thin
	(iSample-1)*thin
}
attr(iSample2time,"ex") <- function(){
	iSample <- 1:13
	structure( iSample2time(iSample) , names=iSample )
	structure( iSample2time(iSample, thin=10) , names=iSample )
	structure( iSample2time(iSample, thin=4) , names=iSample )
	structure( iSample2time(iSample, thin=8) , names=iSample )
}

time2iSample <- function(
	### Convert time to sample number
	time	##<< numeric: indices of the samples in recorded states
	,thin=1	##<< integer scalar: thinning interval
	,match=c(	##<< mode of matching times between samples
		##describe<<
		round="round"		##<< closest sample		
		,floor="floor"		##<< sample before time
		,ceiling="ceiling"	##<< sample after time
		,none="none"		##<< returns a fractional, use if you sure that time corresponds to actual sample indices
		)##end<<
){
	##seealso<<   
	## \code{\link{iSample2time}}
	## \code{\link{getNGen.twDEMC}}
	## ,\code{\link{subChains.twDEMC}}
	## ,\code{\link{twDEMCBlockInt}}
	x <- (time/thin)+1
	match=match.arg(match)
	switch(match
		,round = round(x)
		,floor = floor(x)
		,ceiling = ceiling(x)
		,none = x
	)
	### index of the neares sample
}
attr(time2iSample,"ex") <- function(){
	time = 0:20
	structure( time2iSample(time), names=time)
	structure( iSample2time(1:6, thin=4) , names=1:6 )	# to show the times corresponding to sample
	structure( time2iSample(time,thin=4,match="floor"), names=time)
	structure( time2iSample(time,thin=4,match="ceiling"), names=time)
	structure( time2iSample(time,thin=4), names=time)	# round may vary for times exactly between two samples
	structure( time2iSample(time,thin=4,match="none"), names=time)	# here we do not get indices
}



R.methodsS3::setMethodS3("getNGen","twDEMC", function( 
		### Extract the number of completed generations in res
		res	##<< object of class twDEMC
		,... 
	){
		##details<< 
		## the number of generations corresponds to the steps after time 0 (sample 1).
		## Hence the sample of size 2 and thinning 1 describes one generation (one step forward).
		## A sample of size 2 of thinning 5 encompasses times 0 and 5, i.e. 5 generations.
		## see 
		
		# getNGen.twDEMC
		##seealso<<   
		## \code{\link{iSample2time}}
		## \code{\link{time2iSample}}
		## \code{\link{getNPops.twDEMC}}
		## \code{\link{getNSamples.twDEMC}}
		## \code{\link{getNChains.twDEMC}}
		## \code{\link{getNChainsPop.twDEMC}}
		## \code{\link{getNParms.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCBlockInt}}
		(nrow(res$logDen)-1)*res$thin
		### integer, number of completed generations
	})
attr( getNGen.twDEMC, "ex") <- function(){
	data(twdemcEx1)
	getNGen(twdemcEx1)
	getNSamples(twdemcEx1)
	twdemcEx1$thin
	getNPops(twdemcEx1)
	getNChains(twdemcEx1)
	getNChainsPop(twdemcEx1)
	getNParms(twdemcEx1)
	getNBlocks(twdemcEx1)
}

R.methodsS3::setMethodS3("getNBlocks","twDEMC", function( 
		### Extracts the number of blocks 
		x	##<< object of class twDEMC
		,... 
	){
		# getNPops.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## ,\code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCBlockInt}}
		ncol(x$logDen)
		### integer, number of blocks in twDEMC
	})


R.methodsS3::setMethodS3("getNPops","twDEMC", function( 
		### Extracts the number of populations
		res	##<< object of class twDEMC
		,... 
	){
		# getNPops.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## ,\code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCBlockInt}}
		if( 0==length(res$nPop) ){
			warning("getNPops.twDEMC: mssing $nPop, returning 1")
			1
		}else{
			res$nPop
		}
		#ncol(res$temp)
		### integer, number of populations in twDEMC
	})

R.methodsS3::setMethodS3("getNChains","twDEMC", function( 
		### Extracts the number of chains
		res	##<< object of class twDEMC
		,... 
	){
		# getNChains.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCBlockInt}}
		dim(res$parms)[3]
		### integer, number of chains in twDEMC
	})

R.methodsS3::setMethodS3("getNChainsPop","twDEMC", function( 
		### Extracts the number of chains per population
		res	##<< object of class twDEMC
		,... 
	){
		# getNChainsPop.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCBlockInt}}
		dim(res$parms)[3]/res$nPop
		### integer, number of chains per population in twDEMC
	})

R.methodsS3::setMethodS3("getNParms","twDEMC", function( 
		### Extracts the number of parameters, i.e. the dimension of the parameter vector
		res	##<< object of class twDEMC
		,... 
	){
		# getNParms.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCBlockInt}}
		ncol(res$parms)
		### integer, number of parameters in twDEMC
	})

R.methodsS3::setMethodS3("getNSamples","twDEMC", function( 
		### Extracts the number of s
		res	##<< object of class twDEMC
		,... 
	){
		##details<< There is only one sample per thinning interval of length \code{res$thin}.
		# getNSamples.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCBlockInt}}
		nrow(res$parms)
		### integer, number of samples in twDEMC
	})

R.methodsS3::setMethodS3("getCurrentTemp","twDEMC", function( 
		### Get the Temperature, i.e. cost reduction factor of the last sample
		x	##<< object of class twDEMCPops
		,... 
	){
		# getCurrentTemp.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		## \link{twDEMC}
		##details<< There is only one sample per thinning interval of length \code{x$thin}.
		x$temp[ nrow(x$temp), ]
		### numeric vector: Temperature for each result component for the last sample
	})




#R.methodsS3::setMethodS3("as.mcmc.list","twDEMC", function( 
as.mcmc.list.twDEMC <- function( 
	### Converts list of type twDEMC (result of \code{\link{twDEMC}}) to coda's \code{mcmc.list}. 
	x,				##<< the output of \code{\link{twDEMC}}) run
	maxLength=NULL, ##<< maximum length of the chains, thin will be increased accordingly
	thin=1, 		##<< thinning interval
	start=1,		##<< starting generation (for neglecting the burnin)
	...
){
	# as.mcmc.list.twDEMC
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	origThin <- x$thin
	if( !is.null(maxLength) ){
		l <- (dim(x$parms)[2] - start%/%origThin)*origThin
		mlThin <- ceiling(l/maxLength)
		#must be a multiple of origThin
		mlThin = origThin * ceiling(mlThin/origThin) 
	}else mlThin = 1
	thin = max( origThin, thin, mlThin )
	nChain = dim(x$parms)[3]
	tmp.l <- lapply( 1:nChain, function(i){
			window(	mcmc( (adrop(x$parms[,,i,drop=FALSE],3)), thin=origThin ), start=start, thin=thin )
		})	
	mcmc.list( tmp.l )
	### a \code{\link[coda]{mcmc.list}}
}
#)
#mtrace(as.mcmc.list.twDEMC)


R.methodsS3::setMethodS3("subset","twDEMC", function( 
	### Condenses an twDEMC result object to the cases boKeep.
	x 			##<< twDEMC object
	,boKeep		##<< either logical vector or numeric vector of indices of cases to keep
	,...
){
	# subset.twDEMC 
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	if( is.logical(boKeep) )
		boKeep <- rep(boKeep,length.out=nrow(x$logDen) )
	x$parms <- x$parms[boKeep,,, drop=FALSE] 
	x$logDen <- x$logDen[boKeep,,, drop=FALSE] 
	x$resLogDen <- x$resLogDen[boKeep,,, drop=FALSE] 
	x$pAccept <- x$pAccept[boKeep,,, drop=FALSE]
	x$temp <- x$temp[boKeep, ,drop=FALSE]
	##details<<
	## components \code{thin,Y,nGenBurnin} are kept, but may be meaningless after subsetting.
	x
	### list of class twDEMC with subset of cases in parsm, logDen, pAccept, and temp
})
#mtrace(subset.twDEMC)

R.methodsS3::setMethodS3("thin","twDEMC", function( 
	### Reduces the rows of an twDEMC object (list returned by \code{\link{twDEMCBlockInt}}) to correspond to a thinning of \code{newThin}.
	x, ##<< the twDEMC list to thin 
	newThin=x$thin, ##<< finite numeric scalar: the target thinning factor, must be positive multiple of x$thin 
	start=0,  
		### numeric scalar: the start time of the chain. 
		### Note that time starts from zero.
		### If a vector or matrix is supplied (e.g. nGenBurnin) then the maximum is used
	end=NA,   
		### numeric scalar: the maximum end time of the chains. 
		### Note that time starts from zero.
		### If a vector or matrix is supplied (e.g. nGenBurnin) then the maximum is used
	...
	, doKeepBatchCall=FALSE	##<< wheter to retain the batch call attribute of x

){
	# thin.twDEMC
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
		
	# with the thinned list having mZ rows, this corresponds to (mZ-1)*thin Metropolis steps + 1 row for the initial state
	if( (newThin < x$thin) | (newThin %% x$thin) )
		stop(paste("thin.twDEMC: increased thin must be a positive multiple of former thin",x$thin))
	start <- max(start)	#
	if( start < 0)
		stop(paste("thin.twDEMC: argument start must be at least 0 but was ",start))
	nS <- getNSamples(x)
	maxSampleTime <- iSample2time(nS, thin=x$thin)
	suppressWarnings( end<-min(end) )
	if( is.null(end) || !is.finite(end) || end>maxSampleTime) end=maxSampleTime 	
	if( end < 1)
		stop(paste("thin.twDEMC: argument end must be at least 1 (one generation from 0 to 1) but was",end))
	#thin own past: keep first line and every occurenc of multiple thin
	nGen <- getNGen(x)
	thinFac <- newThin %/% x$thin
	startT <- ceiling( start / x$thin ) * x$thin		# adjust start time so that it coincides with next start of next thinning interval
	endT <- startT + floor( (end-startT) / newThin) * newThin 			# adjust end time so that it coincides with beginning of thinning interval of end
	iStartEnd <- time2iSample( c(startT,endT), thin=x$thin, match="none" )
	iKeep <- seq(iStartEnd[1],iStartEnd[2],by=thinFac)
	res <- subset.twDEMC( x, iKeep )
	res$thin <- newThin
	#time2iSample(70,5)
	if( !is.null(x$nGenBurnin) ) res$nGenBurnin <- pmax(0,x$nGenBurnin-startT) 
	if(!doKeepBatchCall) attr(res,"batchCall") <- NULL
	res
})
#mtrace(thin.twDEMC)
attr(thin.twDEMC,"ex") <- function(){
	data(twdemcEx1)
	x <- twdemcEx1
	c( nGen=getNGen(twdemcEx1), thin=twdemcEx1$thin, nSample=getNSamples(twdemcEx1) )

    .nGenBurnin <- max(getNGen(twdemcEx1))-70
	thinned <- thin(twdemcEx1, start=.nGenBurnin)	# removing burnin period
	c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin )	#18 sample describing 68 generations

	thinned <- thin(twdemcEx1, start=.nGenBurnin, newThin=8)	
	c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin )	#9 samples describing 64 generations
}
#tw



R.methodsS3::setMethodS3("combinePops","twDEMC", function( 
	### Combine several populations to one big population consiting of more chains.
	x		##<< first twDEMC object
	, ...	##<< more twDEMC objects
	, doKeepBatchCall=FALSE	##<< wheter to retain the batch call attribute of x
){
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	##details<< 
	## All arguments must be of type twDEMC and of same thinning interval and length.
	## Temperature is taken from the first population - it assumes that all populations have the same Temperature in their steps
	res <- x
	.pops <- c( list(res), list(...))
	if( !all(sapply(.pops, is, "twDEMC")) )
		stop("combinePops.twDEMC: all arguments must be of class twDEMC")
	.d <- sapply( .pops,function(x)dim(x$parms))
	if( sum(.d-.d[,1])>0 )
		stop("combinePops.twDEMC: all arguments must have the same dimension")
	if( sum(sapply(.pops,function(x)x$thin)-res$thin)>0 )
		stop("combinePops.twDEMC: all arguments must have the same thinning")
	res$parms <- abind( lapply(.pops, function(psr){psr$parms}) )
	if( 0 <length(x$Y) ){
		res$Y <- abind( lapply(.pops, function(psr){psr$Y}) )
		dimnames(res$Y) <- dimnames(x$Y)
	}
	res$logDen <- abind( lapply(.pops, function(psr){psr$logDen}) )
	res$logDenComp <- abind( lapply(.pops, function(psr){psr$logDenComp}) )
	res$nGenBurnin <- sapply(.pops, "[[", "nGenBurnin") #TODO
	if( !is.numeric(res$nGenBurnin) )
		res$nGenBurnin <- NULL
	res$pAccept <- abind( lapply(.pops, function(psr){psr$pAccept}) )
	res$thin <- x$thin
	#res$temp <- abind( lapply(.pops, function(psr){psr$temp}), along=2 )
	res$temp <- x$temp
	for( i in c("parms","logDen","logDenComp","pAccept","temp") )
		dimnames(res[[i]]) <- dimnames(x[[i]])
	if( !doKeepBatchCall ) attr(res,"batchCall") <- NULL	#keep the 
	res
})
#mtrace(combinePops.twDEMC)


R.methodsS3::setMethodS3("stackChains","array", function( 
	### Combine MarkovChains of result population of twDEMC. 
	Zinit,
	...
){
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	nChain <- dim(Zinit)[3]
	resChains <- lapply( 1:nChain, function(i){ Zinit[,,i,drop=FALSE] }) 
	res <- adrop(abind( resChains, along=1 ),3)
	dimnames(res) <- dimnames(Zinit)[1:2]
	res
	### Numeric matrix (nStep x nParm). steps=c(steps_chain1, steps_chain2_, ... )
})
#mtrace(stackChains.array)
#stackChains(Zinit)

#(tres <- twUtestF(combinePops,"test.stackChains"))
R.methodsS3::setMethodS3("stackChains","twDEMC", function( 
	### Combine MarkovChains of a twDEMC to a matrix. 
	x
	,omitBurnin=FALSE	##<< if TRUE, then burnin of each chain is omitted before stacking
	,...
){
	#stackChains.twDEMC
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	nPop <- getNPops(x)
	start <- if( omitBurnin && !is.null(x$nGenBurnin) ) ceiling(x$nGenBurnin)  else rep(0,nPop)  
	if( length(start)==1 ) start=rep(start,nPop)
	if( length(start) != nPop) stop("stackChains.twDEMC: burnin must be of length of number of populations.")
	startChain <- rep(start, each=getNChainsPop(x) )
	nChain <- getNChains(x)
	res <- cbind( logDen = adrop(abind( lapply( 1:nChain, function(i){ 	
					if( startChain[i] == 0) x$logDen[,,i ,drop=FALSE] else x$logDen[-(1:(startChain[i]%/%x$thin)),,i ,drop=FALSE]		}
		), along=1),3 )
		#,parms = stackChains.array(x$parms)
		,parms = adrop(abind( lapply( 1:nChain, function(i){	
							if( startChain[i] == 0) x$parms[,,i ,drop=FALSE] else x$parms[-(1:(startChain[i]%/%x$thin)),,i ,drop=FALSE]	
						}), along=1 ),3)
	)
	attr(res,"nBlock") = getNBlocks(x)
	res
	### Matrix with first attributes$nBlock columns the logDensity logDen and the remaining columns the variables.
})

R.methodsS3::setMethodS3("stackChains","mcmc.list", function( 
	### Combine list of MarkovChains to one big matrix.
	mcmcList,
	...
){
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	abind(mcmcList, along=1)
})


mcmcListApply <- function(
	### Applies function to each chain of \code{mcmc.list}.
	mcmcList,
	fun,
	...
){
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	res <- lapply( mcmcList, fun, ...)
	attributes(res) <- attributes(mcmcList)
	res
}



#mtrace(plotThinned.mcmc.list)


popMeansTwDEMC <- function( 
	### Calculating population means across chains within each population, and smooth time series.
	x				##<< a matrix with columns chains
	,nPop			##<< number of populations 
	,kSmooth=NULL	##<< weights to the filter function, or just number of points
){
	#popMeansTwDEMC
	##seealso<<   
	## \code{\link{popApplyTwDEMC}}
	## \code{\link{subChains.twDEMC}}
		
	boMat <- is.matrix(x)
	if( !boMat) x <- matrix(x,nrow=1)
	p <- list( nrow=nrow(x), ncol=ncol(x) )
	nChainsP <- p$ncol %/% nPop
	xPop <- matrix( double(p$nrow*nPop), ncol=nPop, dimnames=c( dimnames(x)[1], list(pops=NULL)) )
	for( iPop0 in (0:(nPop-1)) ){
		iChains <-  iPop0*nChainsP + (1:nChainsP)
		xPop[,iPop0+1] <- rowMeans(x[,iChains,drop=FALSE])		
	}
	if( is.numeric(kSmooth) ){
		if( length(kSmooth)==1) kSmooth=rep(1/kSmooth,kSmooth)
		xPop <- filter(xPop,kSmooth)
	}
	### matrix of dim( nrow(x) x nPop ) with rowMeans across chains of one population.
	if( boMat) xPop else as.vector(xPop)
}
attr(popMeansTwDEMC,"ex") <- function(){
	data(twdemcEx1)
	ex1c <- concatPops(twdemcEx1)
	# mean logDen for each case, i.e. step, by population
	#mtrace(popMeansTwDEMC)
	res1 <- popMeansTwDEMC( ex1c$logDen[,1,], nPop=getNPops(ex1c) )
	matplot(res1)
	# shifting mean across 4 cases
	res2 <- popMeansTwDEMC( ex1c$logDen[,1,], nPop=getNPops(ex1c), kSmooth=5 )
	matplot(res2, type="l", add=TRUE)
}

popApplyTwDEMC <- function( 
	### Applying a function across all chains of one population for each case.
	x				##<< a matrix with columns chains or array with last dimension chain
	,nPop			##<< number of populations
	,FUN			##<< function to apply to population submatrix
	,...			##<< further arguemtns to FUN
){
	#popApplyTwDEMC
	##seealso<<   
	## \code{\link{popMeansTwDEMC}}
	## \code{\link{subChains.twDEMC}}
	ldim <-length(dim(x)) 
	nChainsP <- dim(x)[ldim] %/% nPop
	resl <- 
		if( ldim==2 ) 
	 		lapply( (0:(nPop-1)), function(iPop0){
				iChains <-  iPop0*nChainsP + (1:nChainsP)
				FUN(x[,iChains ,drop=FALSE],...) 		
			})
		else if( ldim==3 )
			#iPop0=0
			lapply( (0:(nPop-1)), function(iPop0){
					iChains <-  iPop0*nChainsP + (1:nChainsP)
					FUN(x[,,iChains ,drop=FALSE],...) 		
				})
		else stop("dimension of x must be 2 or 3")
	abind(resl, rev.along=0)
	#resl
	### array with last dimenstion correponding to population
}
attr(popApplyTwDEMC,"ex") <- function(){
	data(twdemcEx1)
	ex1c <- concatPops(twdemcEx1)
	# mean logDen for each case, i.e. step, by population
	nPop=getNPops(ex1c)
	#mtrace(popApplyTwDEMC)
	# applied to a matrix: pops in columns
	popApplyTwDEMC( ex1c$logDen[,1,], nPop=nPop, apply, 1, mean )
	# applied to a 3d array: pops in 3rd dimension
	tmp <- popApplyTwDEMC( ex1c$logDen, nPop=nPop, apply, 1:2, mean )
	str(tmp)
	
	# stack logDen of block 1 for each population
	tmp <- popApplyTwDEMC( ex1c$logDen[,1,], nPop=nPop, as.vector )
	str(tmp)
	#stack param columns by population
	str(ex1c$parms)	#26 cases, 2 parameters, 2*4=8 chains
	tmp <- popApplyTwDEMC( ex1c$parms, nPop=nPop, function(x){ 
                abind::abind(twMisc::twListArrDim(x),along=1) })
	all.equal( c(nrow(ex1c$parms)*4, 2,2), dim(tmp))
}

R.methodsS3::setMethodS3("stackChainsPop","twDEMC", function( 
		### Combine MarkovChains of each population of a twDEMC. 
		x
		,...
	){
		#stackChainsPop.twDEMC
		##seealso<<   
		## \code{\link{stackChains.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		# stack logDen for each population
		nPop = getNPops(x)
		logDen <- popApplyTwDEMC( x$logDen, nPop=nPop, function(x){ abind(twMisc::twListArrDim(x),along=1) })
		#logDen <- popApplyTwDEMC( x$logDen, nPop=nPop, as.vector )
		parms <- popApplyTwDEMC( x$parms, nPop=nPop, function(x){ abind(twMisc::twListArrDim(x),along=1) })		
		res <- abind(logDen,parms,along=2)
		attr(res,"nBlock") <- ncol(logDen)
		res
		### Array with first column the logDensity logDen and the remaining columns the variables
		### , rows are steps, third dimension is the population 
		### (but see argument \code{varInRows}
	})
attr(stackChainsPop.twDEMC,"ex") <- function(){
	data(twdemcEx1)
	ex1c <- concatPops(twdemcEx1)
	res <- stackChainsPop(ex1c)
	str(res)
	(tmp1 <- head(res[,,1]))
}

R.methodsS3::setMethodS3("thinN","mcmc.list", function( 
	### Thin x so that each chain consists of about nThinnedRecords
	x		##<< object to be thinned and plotted (mcmc.list)
	,nThinnedRecords=100	##<< number of records in thinned object
	,...
){
	##details<< 
	## Plotting mcmc.list with large n (>1000) takes long and involves much overplotting.
	## The same visual impression can be achieved using a thinned chain.
	window(x, thin= thin(x)*max(1,(nrow(x[[1]]) %/% nThinnedRecords)))
	### Thinned x
})



.tmp.f <- function(){
	#not adapted to block sampling yet
	R.methodsS3::setMethodS3("replacePops","twDEMC", function( 
			### replaces several populations by other ones
			x
			,xNew
			,iPops 
			,... 
		){
			##seealso<<   
			## \code{\link{subChains.twDEMC}}
			## ,\code{\link{twDEMCBlockInt}}
			
			nPop=getNPops(x)
			nPopNew = getNPops(xNew)
			nChains = ncol(x$logDen)
			nChainPop = nChains %/% nPop
			if( length(iPops) != nPopNew) stop("replacePops.twDEMC: replacing population has number of populations that differs from length of iPop.")
			nChainPopNew = ncol(xNew$logDen) %/% nPopNew
			if( nChainPop != nChainPopNew) stop("replacePops.twDEMC: both populations must have the same number of chains per populations.")
			xN <- if( x$thin != xNew$thin ) thin(xNew,newThin=x$thin) else xNew
			if( !identical( dim(x)[1:2], dim(xN)[1:2]) ) stop("replacePops.twDEMC: both populations must have the same dimensions of variables and cases.")
			
			res <- x
			iChains <- matrix(1:nChains, nrow=nChainPop)[,iPops]
			#now no differenc ein T by populations res$temp[,iPops] <- xN$temp  
			res$parms[,,iChains] <- xN$parms
			res$logDenComp[,,iChains] <- xN$logDenComp
			res$logDen[,iChains] <- xN$logDen 
			res$pAccept[,iChains] <- xN$pAccept
			if( 0 < length(x$Y))
				res$Y[,,iChains] <- xN$Y
			res
			### a list of class twDEMC (see \code{\link{twDEMCBlockInt}})
		})
#mtrace(subChains.twDEMC)
}



