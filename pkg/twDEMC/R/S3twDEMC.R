#	(tres <- twUtestF("S3twDEMC")) #thin, subset, subChains, combinePops

setMethodS3("subChains","twDEMC", function( 
		### Condenses an twDEMC List to the chains iChains e.g \code{1:4}.
		x, iChains=NULL, iPops=NULL, 
		nPop=NULL, ##<< number of populations in x, if not specified then taken from ncol(x$temp)
		... 
		, doKeepBatchCall=FALSE	##<< wheter to retain the batch call attribute of x
){
		##seealso<<   
		## \code{\link{twDEMCInt}}
		
		##details<< 
		## There are several methods access properties a result of an \code{\link{twDEMCInt}}, i.e.
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
		## There are several methods to transform or subset the results of an \code{\link{twDEMCInt}} run. \itemize{
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
		## \item{ plotting a subset of the chains and cases: \code{\link{plotThinned.mcmc.list}}  } 
		## \item{ transforming parameters \code{\link{transOrigPopt.mcmc.list}}  } 
		##}

		# To help re-initializing the arguments to fLogDen \itemize{
		# \item{ transforming parameters \code{\link{subsetArgsFLogDen}}  }}
		#
		
		##details<< 
		## Alternatively to specification of iChains, one can specify a vector of populations and the total number of populations.
		if( is.null(nPop) ) nPop=ncol(x$temp)
		if( is.null(nPop) ) nPop=1		#for backware compatibility of temp
		nChainsPop = ncol(x$rLogDen) %/% nPop
		res <- x
		if( !is.null(iPops)){
			iChains = unlist( lapply( iPops, function(iPop){ (iPop-1)*nChainsPop + (1:nChainsPop) }))
			res$temp <- x$temp[,iPops, drop=FALSE]	#by nPops
			if( !is.null(x$nGenBurnin) ) res$nGenBurnin <- x$nGenBurnin[iPops]
		}else{
			# no populatios given: reduce to one population
			res$temp <- try( matrix( x$temp[,1], nrow=nrow(x$temp), ncol=1 ),silent=TRUE)
			if( inherits(res$temp,"try-error")) res$temp=matrix( 1, nrow=nrow(x$rLogDen), ncol=1) #for backware compatibility of temp
			if( !is.null(x$nGenBurnin) ) res$nGenBurnin <- max(x$nGenBurnin)
		}
		res$parms <- x$parms[,,iChains, drop=FALSE]
		res$logDenComp <- x$logDenComp[,,iChains, drop=FALSE]
		res$rLogDen <- x$rLogDen[,iChains, drop=FALSE] 
		res$pAccept <- x$pAccept[,iChains, drop=FALSE]
		res$thin <- x$thin
		if( 0 < length(x$Y))
			res$Y <- x$Y[,,iChains, drop=FALSE]
		res$thin <- x$thin
		if(!doKeepBatchCall) attr(res,"batchCall") <- NULL
		res
		### a list of class twDEMC (see \code{\link{twDEMCInt}})
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
	## ,\code{\link{twDEMCInt}}
	
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
	## ,\code{\link{twDEMCInt}}
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
attr(iSample2time,"ex") <- function(){
	time = 0:20
	structure( time2iSample(time), names=time)
	structure( iSample2time(1:6, thin=4) , names=1:6 )	# to show the times corresponding to sample
	structure( time2iSample(time,thin=4,match="floor"), names=time)
	structure( time2iSample(time,thin=4,match="ceiling"), names=time)
	structure( time2iSample(time,thin=4), names=time)	# round may vary for times exactly between two samples
	structure( time2iSample(time,thin=4,match="none"), names=time)	# here we do not get indices
}



setMethodS3("getNGen","twDEMC", function( 
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
		## ,\code{\link{twDEMCInt}}
		(nrow(res$rLogDen)-1)*res$thin
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
}

setMethodS3("getNPops","twDEMC", function( 
		### Extracts the number of populations
		res	##<< object of class twDEMC
		,... 
	){
		# getNPops.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## ,\code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCInt}}
		ncol(res$temp)
		### integer, number of populations in twDEMC
	})

setMethodS3("getNChains","twDEMC", function( 
		### Extracts the number of chains
		res	##<< object of class twDEMC
		,... 
	){
		# getNChains.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCInt}}
		ncol(res$rLogDen)
		### integer, number of chains in twDEMC
	})

setMethodS3("getNChainsPop","twDEMC", function( 
		### Extracts the number of chains per population
		res	##<< object of class twDEMC
		,... 
	){
		# getNChainsPop.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCInt}}
		ncol(res$rLogDen)/ncol(res$temp)
		### integer, number of chains per population in twDEMC
	})

setMethodS3("getNParms","twDEMC", function( 
		### Extracts the number of parameters, i.e. the dimension of the parameter vector
		res	##<< object of class twDEMC
		,... 
	){
		# getNParms.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCInt}}
		nrow(res$parms)
		### integer, number of parameters in twDEMC
	})

setMethodS3("getNSamples","twDEMC", function( 
		### Extracts the number of s
		res	##<< object of class twDEMC
		,... 
	){
		##details<< There is only one sample per thinning interval of length \code{res$thin}.
		# getNSamples.twDEMC
		##seealso<<   
		## \code{\link{getNGen.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCInt}}
		ncol(res$parms)
		### integer, number of samples in twDEMC
	})



setMethodS3("replacePops","twDEMC", function( 
		### replaces several populations by other ones
		x
		,xNew
		,iPops 
		,... 
	){
		##seealso<<   
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCInt}}
		
		nPop=ncol(x$temp)
		nPopNew = ncol(xNew$temp)
		nChains = ncol(x$rLogDen)
		nChainsPop = nChains %/% nPop
		if( length(iPops) != nPopNew) stop("replacePops.twDEMC: replacing population has number of populations that differs from length of iPop.")
		nChainsPopNew = ncol(xNew$rLogDen) %/% nPopNew
		if( nChainsPop != nChainsPopNew) stop("replacePops.twDEMC: both populations must have the same number of chains per populations.")
		xN <- if( x$thin != xNew$thin ) thin(xNew,newThin=x$thin) else xNew
		if( !identical( dim(x)[1:2], dim(xN)[1:2]) ) stop("replacePops.twDEMC: both populations must have the same dimensions of variables and cases.")

		res <- x
		iChains <- matrix(1:nChains, nrow=nChainsPop)[,iPops]
		res$temp[,iPops] <- xN$temp
		res$parms[,,iChains] <- xN$parms
		res$logDenComp[,,iChains] <- xN$logDenComp
		res$rLogDen[,iChains] <- xN$rLogDen 
		res$pAccept[,iChains] <- xN$pAccept
		if( 0 < length(x$Y))
			res$Y[,,iChains] <- xN$Y
		res
		### a list of class twDEMC (see \code{\link{twDEMCInt}})
	})
#mtrace(subChains.twDEMC)


#setMethodS3("as.mcmc.list","twDEMC", function( 
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


setMethodS3("subset","twDEMC", function( 
	### Condenses an twDEMC result object to the cases boKeep.
	x 			##<< twDEMC object
	,boKeep		##<< either logical vector or numeric vector of indices of cases to keep
	,...
){
	# subset.twDEMC 
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	if( is.logical(boKeep) )
		boKeep <- rep(boKeep,length.out=nrow(x$rLogDen) )
	x$parms <- x$parms[,boKeep,, drop=FALSE] 
	x$rLogDen <- x$rLogDen[boKeep,, drop=FALSE] 
	x$logDenComp <- x$logDenComp[,boKeep,, drop=FALSE] 
	x$pAccept <- x$pAccept[boKeep,, drop=FALSE]
	if( !is.null(x$temp))
		if( is.null(ncol(x$temp)) )
			x$temp <- x$temp[boKeep, drop=FALSE]
		else
			x$temp <- x$temp[boKeep,, drop=FALSE]
	##details<<
	## components \code{thin,Y,nGenBurnin} are kept, but may be meaningless after subsetting.
	x
	### list of class twDEMC with subset of cases in parsm, rLogDen, pAccept, and temp
})
#mtrace(subset.twDEMC)

setMethodS3("thin","twDEMC", function( 
	### Reduces the rows of an twDEMC object (list returned by \code{\link{twDEMCInt}}) to correspond to a thinning of \code{newThin}.
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
	startT <- ceiling( start / newThin ) * newThin		# adjust start time so that it coincides with next start of next thinning interval
	endT <- floor( end / newThin) * newThin 			# adjust end time so that it coincides with beginning of thinning interval of end
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
	c( nGen=getNGen(twdemcEx1), thin=twdemcEx1$thin, nSample=getNSamples(twdemcEx1), nGenBurnin=twdemcEx1$nGenBurnin )

	thinned <- thin(twdemcEx1, start=twdemcEx1$nGenBurnin)	# removing burnin period
	c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin )	#15 sample describing 70 generations

	thinned <- thin(twdemcEx1, start=twdemcEx1$nGenBurnin, newThin=10)	
	c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin )	#8 samples describing 70 generations
}
#tw



setMethodS3("combinePops","twDEMC", function( 
	### Combine several populations to one big population.
	x		##<< first twDEMC object
	, ...	##<< more twDEMC objects
	, doKeepBatchCall=FALSE	##<< wheter to retain the batch call attribute of x
){
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	##details<< 
	## All arguments must be of type twDEMC and of same thinning interval and length.
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
	res$rLogDen <- abind( lapply(.pops, function(psr){psr$rLogDen}) )
	res$logDenComp <- abind( lapply(.pops, function(psr){psr$logDenComp}) )
	res$nGenBurnin <- sapply(.pops, "[[", "nGenBurnin") 
	res$pAccept <- abind( lapply(.pops, function(psr){psr$pAccept}) )
	res$thin <- x$thin
	res$temp <- abind( lapply(.pops, function(psr){psr$temp}), along=2 )
	for( i in c("parms","rLogDen","logDenComp","pAccept","temp") )
		dimnames(res[[i]]) <- dimnames(x[[i]])
	if( !doKeepBatchCall ) attr(res,"batchCall") <- NULL	#keep the 
	res
})
#mtrace(combinePops.twDEMC)


setMethodS3("stackChains","array", function( 
	### Combine MarkovChains of a initializer of twDEMC. 
	Zinit,
	...
){
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	res <- t(adrop(abind( lapply( 1:dim(Zinit)[3], function(i){ Zinit[,,i,drop=FALSE] }), along=2 ),3))
	res
	### Matrix with columns the variables.
})
#mtrace(stackChains.array)
#stackChains(Zinit)

#(tres <- twUtestF(combinePops,"test.stackChains"))
setMethodS3("stackChains","twDEMC", function( 
	### Combine MarkovChains of a twDEMC to a matrix. 
	x
	,omitBurnin=FALSE	##<< if TRUE, then burnin of each chain is omitted before stacking
	,...
){
	#stackChains.twDEMC
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	nPop <- getNPops(x)
	start <- if( omitBurnin & !is.null(x$nGenBurnin) ) ceiling(x$nGenBurnin)  else rep(0,nPop)  
	if( length(start)==1 ) start=rep(start,nPop)
	if( length(start) != nPop) stop("stackChains.twDEMC: burnin must be of length of number of populations.")
	startChain <- rep(start, each=getNChainsPop(x) )
	nChain <- getNChains(x)
	cbind( rLogDen = abind( lapply( 1:nChain, function(i){ 	if( startChain[i] == 0) x$rLogDen[,i] else x$rLogDen[-(1:(startChain[i]%/%x$thin)),i]		}), along=1 )
		#,parms = stackChains.array(x$parms)
		,parms = t(adrop(abind( lapply( 1:nChain, function(i){	if( startChain[i] == 0) x$parms[,,i,drop=FALSE] else x$parms[,-(1:(startChain[i]%/%x$thin)),i,drop=FALSE]	}), along=2 ),3))
	)
	### Matrix with first column the logDensity rLogDen and the remaining columns the variables.
})

setMethodS3("stackChains","mcmc.list", function( 
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
	,nPops			##<< number of populations 
	,kSmooth=NULL	##<< weights to the filter function, or just number of points
){
	#popMeansTwDEMC
	##seealso<<   
	## \code{\link{popApplyTwDEMC}}
	## \code{\link{subChains.twDEMC}}
		
	boMat <- is.matrix(x)
	if( !boMat) x <- matrix(x,nrow=1)
	p <- list( nrow=nrow(x), ncol=ncol(x) )
	nChainsP <- p$ncol %/% nPops
	xPop <- matrix( double(p$nrow*nPops), ncol=nPops ) 
	for( iPop0 in (0:(nPops-1)) ){
		iChains <-  iPop0*nChainsP + (1:nChainsP)
		xPop[,iPop0+1] <- rowMeans(x[,iChains,drop=FALSE])		
	}
	if( is.numeric(kSmooth) ){
		if( length(kSmooth)==1) kSmooth=rep(1/kSmooth,kSmooth)
		xPop <- filter(xPop,kSmooth)
	}
	### matrix of dim( nrow(x) x nPops ) with rowMeans across chains of one population.
	if( boMat) xPop else as.vector(xPop)
}
attr(popMeansTwDEMC,"ex") <- function(){
	data(twdemcEx1)
	# mean rLogDen for each case, i.e. step, by population
	res1 <- popMeansTwDEMC( twdemcEx1$rLogDen, nPops=ncol(twdemcEx1$temp) )
	matplot(res1)
	# shifting mean across 4 cases
	res2 <- popMeansTwDEMC( twdemcEx1$rLogDen, nPops=ncol(twdemcEx1$temp), kSmooth=5 )
	matplot(res2, type="l", add=TRUE)
}

popApplyTwDEMC <- function( 
	### Applying a function across all chains of one population for each case.
	x				##<< a matrix with columns chains or array with last dimension chain
	,nPops			##<< number of populations
	,FUN			##<< function to apply to population submatrix
	,...			##<< further arguemtns to FUN
){
	#popApplyTwDEMC
	##seealso<<   
	## \code{\link{popMeansTwDEMC}}
	## \code{\link{subChains.twDEMC}}
	ldim <-length(dim(x)) 
	nChainsP <- dim(x)[ldim] %/% nPops
	resl <- 
		if( ldim==2 ) 
	 		lapply( (0:(nPops-1)), function(iPop0){
				iChains <-  iPop0*nChainsP + (1:nChainsP)
				FUN(x[,iChains,drop=FALSE],...) 		
			})
		else if( ldim==3 ) 
			lapply( (0:(nPops-1)), function(iPop0){
					iChains <-  iPop0*nChainsP + (1:nChainsP)
					FUN(x[,,iChains,drop=FALSE],...) 		
				})
		else stop("dimension of x must be 2 or 3")
	abind(resl, rev.along=0)
	#resl
	### array with last dimenstion correponding to population
}
attr(popApplyTwDEMC,"ex") <- function(){
	data(twdemcEx1)
	# mean rLogDen for each case, i.e. step, by population
	nPops=getNPops(twdemcEx1)
	popApplyTwDEMC( twdemcEx1$rLogDen, nPops=nPops, apply, 1, mean )	
	# stack rLogDen for each population
	popApplyTwDEMC( twdemcEx1$rLogDen, nPops=nPops, as.vector )
	#stack param columns by population
	(tmp <- popApplyTwDEMC( twdemcEx1$parms, nPops=nPops, function(x){ abind(twListArrDim(x),along=2) }))	
}

setMethodS3("stackChainsPop","twDEMC", function( 
		### Combine MarkovChains of each population of a twDEMC. 
		x
		,...
		,varInRows=FALSE	##<< set to TRUE if rows hold variables and columns steps (as in Zinit of twDEMC), defaults to variables in columns
	){
		#stackChainsPop.twDEMC
		##seealso<<   
		## \code{\link{stackChains.twDEMC}}
		## \code{\link{subChains.twDEMC}}
		# stack rLogDen for each population
		nPops = getNPops(x)
		rLogDen <- popApplyTwDEMC( x$rLogDen, nPops=nPops, as.vector )
		#stack param columns by population
		if( varInRows ){
			rLogDen3 <- array( rLogDen, dim=c(1, nrow(rLogDen), ncol(rLogDen) ))
			parms <- popApplyTwDEMC( x$parms, nPops=nPops, function(xp){ abind(twListArrDim(xp),along=2) })	
			res <- structure( abind(rLogDen3,parms,along=1), dimnames=list(parms=c("rLogDen",rownames(x$parms)),step=NULL,pops=NULL))
		}else{
			rLogDen3 <- array( rLogDen, dim=c(nrow(rLogDen),1,ncol(rLogDen)))
			parms <- popApplyTwDEMC( x$parms, nPops=nPops, function(xp){ t(abind(twListArrDim(xp),along=2)) })	
			res <- structure( abind(rLogDen,parms,along=2), dimnames=list(step=NULL,parms=c("rLogDen",rownames(x$parms)),pops=NULL))
		}
		res
		### Array with first column the logDensity rLogDen and the remaining columns the variables
		### , rows are steps, third dimension is the population 
		### (but see argument \code{varInRows}
	})
attr(stackChainsPop.twDEMC,"ex") <- function(){
	data(twdemcEx1)
	res <- stackChainsPop(twdemcEx1)
	str(res)
	(tmp1 <- head(res[,,1]))
	
	res2 <- stackChainsPop(twdemcEx1, varInRows=TRUE)
	str(res2)
	(tmp2 <- head(t(res2[,,1])))
	
	identical( tmp1, tmp2)
}





setMethodS3("thinN","mcmc.list", function( 
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

getAcceptedPos.twDEMCProps <- function(
	### Calculate the accepted state corresponding to the proposal.
	Y	##<< matrix of proposals with row "accepted" and first step (column) initial state, and third dimension chains.
	,acceptedRowName = "accepted"
){
	#j=1
	# set the accepted state of the first entry to TRUE so that a start is given 
	Y[acceptedRowName,1,] <- 1
	tmp.f <- function(j){		#j is the chain, i.e. column
		acc <- which( Y[acceptedRowName,,j,drop=TRUE] != 0)
		accl <- diff(c(acc, ncol(Y)+1))
		rep(acc,accl)
	}	
	acceptedPos <- sapply(1:dim(Y)[3], tmp.f) 
	#acceptPos for an accepted proposal is the proposal step itselv, we need the previous one
	if( is.matrix(acceptedPos) )
		prevAcceptedPos <- acceptedPos[c(1,1:(nrow(acceptedPos)-1)),]
	else
		# sapply reduces the case for only 1 step to a vector
		prevAcceptedPos <- matrix(acceptedPos,nrow=1)
	prevAcceptedPos[1,] <- NA
	prevAcceptedPos
	### matrix of steps (steps x chains) whose state was the accepted state when evaluating proposal (row of Y)
	### the first step is NA because there is no previous accepted state
}

getDiffLogDen.twDEMCProps <- function(
	### Extract the Differences in LogDensity between accepted states and proposals.
	Y					##<< matrix of proposals with row "accepted" and first step (column) initial state, rows: results components of fLogDen and third dimension chains.
	,resCols			##<< the rows of Y with result components, either names or positions
	,nLastSteps = 128	##<< number of last steps of Y for which to extract diffs
	,temp=1				##<< numeric matrix (resComp x pop): the temperature applied to the difference of LogDens
	,...				##<< further arguments passed to \code{\link{getAcceptedPos.twDEMCProps}}
){
	#work with positions
	if( is.character(resCols) ) resCols <- match(resCols,rownames(Y))
	# make sure that temp has all resCols components 
	if( 1 < length(temp)) temp <- temp[rownames(Y)[resCols],] else temp <- matrix(temp,ncol=dim(Y)[3],nrow=1)
	nPopsChain <- dim(Y)[3] %/% ncol(temp)
	# constrain to the nLastSteps+1 rows
	YL <- if( nLastSteps+1 >= ncol(Y)) Y else Y[,ncol(Y)+1-((nLastSteps+1):1),,drop=FALSE]
	#test 1d case: YL <- Y[,min(ncol(Y),nLastSteps+1),,drop=FALSE]
	#mtrace(getAcceptedPos.twDEMCProps)
	acceptedPos <- getAcceptedPos.twDEMCProps(YL,...) #index of parameter vector that correponds to the currently accepted for vector at given row
	diffLogDen <- abind(lapply( 1:dim(YL)[3], function(j){ 	adrop((YL[resCols,-1,j,drop=FALSE] - YL[resCols,acceptedPos[-1,j],j,drop=FALSE])/temp[,(j-1)%/%nPopsChain+1],3)}),rev.along=0)
	diffLogDen
	### numeric array ( component x nLastSteps x chain ) of Lp-La, the L
}

replaceNonFiniteDiffLogDens <- function(
	### For each component replace NAs by sample of others and non-Finite values by minimum of others  
	diffLogDen			##<< array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}
	,doConstrainNeg=FALSE	##<< if given, density of accepted jumps (positive) is constrained to 0
){
	d <- t(apply( diffLogDen,1,function(ds){
				boNA <- is.na(ds)
				nNA <- sum(boNA)
				if(0<nNA ){
					dsNonNA <- ds[!boNA]
					ds[boNA] <- sample( dsNonNA, nNA, replace=TRUE)
				}
				boNonFinite <- !is.finite(ds)
				if( 0<sum(boNonFinite)){
					ds[boNonFinite] <- min(ds[!boNonFinite])
				}
				ds
			}))
	dimnames(d) <- dimnames(diffLogDen)
	if( doConstrainNeg )
		d[d>0] <- 0
	d
}

sampleAcceptedFixedTempDiffLogDens <- function(
	### Calculate diffLogDen for fixed temperature components and return subset of accepted ones 
	diffLogDen			##<< array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}
	,TFix=numeric(0)
){
	dTFix <- colSums( diffLogDen[names(TFix),,drop=FALSE]/TFix )
	acceptedFixed <- dTFix>log(runif(ncol(diffLogDen)))
	#pa <- sum(acceptedFixed)/nj
	diffLogDen[,acceptedFixed,drop=FALSE]
}

