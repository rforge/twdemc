

setMethodS3("getNGen","twDEMCPops", function( 
		### Extract the number of completed generations in res
		res	##<< object of class twDEMCPops
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
		## \code{\link{subChains.twDEMCPops}}
		## ,\code{\link{twDEMCInt}}
		#mtrace(getNSamples.twDEMCPops)
		(getNSamples(res)-1)*res$thin
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
		res	##<< object of class twDEMCPops
		,... 
	){
		# getNPops.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## ,\code{\link{subChains.twDEMCPops}}
		## ,\code{\link{twDEMCInt}}
		length(res$pops)
		### integer, number of populations in twDEMCPops
	})

setMethodS3("getNChains","twDEMCPops", function( 
		### Extracts the number of chains
		res	##<< object of class twDEMCPops
		,... 
	){
		# getNChains.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subChains.twDEMCPops}}
		## ,\code{\link{twDEMCInt}}
		length(res$pops) * dim(res$pops[[1]]$parms)[3]
		### integer, number of chains in twDEMCPops
	})

setMethodS3("getNChainsPop","twDEMCPops", function( 
		### Extracts the number of chains per population
		res	##<< object of class twDEMCPops
		,... 
	){
		# getNChainsPop.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subChains.twDEMCPops}}
		## ,\code{\link{twDEMCInt}}
		dim(res$pops[[1]]$parms)[3]
		### integer, number of chains per population in twDEMCPops
	})

setMethodS3("getNParms","twDEMCPops", function( 
		### Extracts the number of parameters, i.e. the dimension of the parameter vector
		res	##<< object of class twDEMCPops
		,... 
	){
		# getNParms.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subChains.twDEMCPops}}
		## ,\code{\link{twDEMCInt}}
		ncol(res$pops[[1]]$parms)
		### integer, number of parameters in twDEMCPops
	})

setMethodS3("getNSamples","twDEMCPops", function( 
		### Extracts the number of samples
		res	##<< object of class twDEMCPops
		,... 
	){
		# getNSamples.twDEMCPops
		##seealso<<   
		## \code{\link{getNGen.twDEMCPops}}
		## \code{\link{subChains.twDEMCPops}}
		## ,\code{\link{twDEMCInt}}
		##details<< There is only one sample per thinning interval of length \code{res$thin}.
		sapply( res$pops, function(pop){ nrow(pop$parms) })
		### integer vector: number of samples in each population of twDEMCPops
	})



setMethodS3("stackPops","twDEMCPops", function( 
	### replaces several populations by other ones
	x
	,... 
	, useThinning=TRUE	##<< if TRUE thinning is used to make populations the same length, if FALSE they are cut to shortest population	
){
	#stackPops.twDEMCPops
	##seealso<<   
	## \code{\link{subChains.twDEMCPops}}
	## ,\code{\link{twDEMCInt}}
	nStepsPop <- getNSamples(x)
	nSteps <- min(nStepsPop)
	if( !all(nStepsPop == nSteps) ){
		if( useThinning)
			x <- thin(x, length.out=nSteps )
		else
			x <- subset(x, 1:nSteps)
	}
	pops <- x$pops
	p1 <- pops[[1]]
	x$pops <- NULL
	x$parms <- structure( abind( lapply(pops,"[[","parms"), along=3), dimnames=dimnames(p1$parms))
	x$temp <- structure( abind( lapply(pops,"[[","temp"), rev.along=0), dimnames=list(steps=NULL,pops=NULL) )
	x$pAccept <- structure( abind( lapply(pops,"[[","pAccept"), along=3), dimnames=dimnames(p1$pAccept))
	x$resLogDen <- structure( abind( lapply(pops,"[[","resLogDen"), along=3), dimnames=dimnames(p1$resLogDen))
	x$logDen <- structure( abind( lapply(pops,"[[","logDen"), along=3), dimnames=dimnames(p1$logDen))
	x$Y <- structure( abind( lapply(pops,"[[","Y"), along=3), dimnames=dimnames(p1$Y))
	class(x) <- c("list","twDEMC")
	x
})

#mtrace(stackPops.twDEMCPops)


#setMethodS3("as.mcmc.list","twDEMCPops", function( 
as.mcmc.list.twDEMCPops <- function( 
	### Converts list of type twDEMCPops (result of \code{\link{twDEMCPops}}) to coda's \code{mcmc.list}. 
	x				##<< the output of \code{\link{twDEMCBlockInt}}) run
	,...
	,useThinning=TRUE	##<< if TRUE thinning is used to make populations the same length, if FALSE they are cut to shortest population	
){
	xStack <- stackPops.twDEMCPops(x, useThinning=useThinning)
	as.mcmc.list.twDEMC(xStack)
}

