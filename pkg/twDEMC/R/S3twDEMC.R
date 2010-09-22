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

		# To help re-initializing the arguments to fLoglik \itemize{
		# \item{ transforming parameters \code{\link{subsetArgsFLogLik}}  }}
		#
		
		##details<< 
		## Alternatively to specification of iChains, one can specify a vector of populations and the total number of populations.
		if( is.null(nPop) ) nPop=ncol(x$temp)
		if( is.null(nPop) ) nPop=1		#for backware compatibility of temp
		nChainsPop = ncol(x$rLogLik) %/% nPop
		res <- x
		if( !is.null(iPops)){
			iChains = unlist( lapply( iPops, function(iPop){ (iPop-1)*nChainsPop + (1:nChainsPop) }))
			res$temp <- x$temp[,iPops, drop=FALSE]	#by nPops
		}else{
			# no populatios given: reduce to one population
			res$temp <- try( matrix( x$temp[,1], nrow=nrow(x$temp), ncol=1 ),silent=TRUE)
			if( inherits(res$temp,"try-error")) res$temp=matrix( 1, nrow=nrow(x$rLogLik), ncol=1) #for backware compatibility of temp
		}
		res$parms <- x$parms[,,iChains, drop=FALSE]
		res$rLogLik <- x$rLogLik[,iChains, drop=FALSE] 
		res$pAccept <- x$pAccept[,iChains, drop=FALSE]
		res$thin <- x$thin
		if( 0 < length(x$resFLogLikX))
			res$resFLogLikX <- x$resFLogLikX[,iChains, drop=FALSE]
		if( 0 < length(x$Y))
			res$Y <- x$Y[,,iChains, drop=FALSE]
		res$thin <- x$thin
		if(!doKeepBatchCall) attr(res,"batchCall") <- NULL
		res
		### a list of class twDEMC (see \code{\link{twDEMCInt}})
	})
#mtrace(subChains.twDEMC)

setMethodS3("calcNGen","twDEMC", function( 
		### Calculates the number of completed generations in res
		res	##<< object of class twDEMC
		,... 
	){
		# calcNGen.twDEMC
		##seealso<<   
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCInt}}
		(nrow(res$rLogLik)-1)*res$thin
		### integer, number of completed generations
	})

setMethodS3("getNPops","twDEMC", function( 
		### Extracts the number of populations
		res	##<< object of class twDEMC
		,... 
	){
		# calcNGen.twDEMC
		##seealso<<   
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCInt}}
		ncol(res$temp)
		### integer, number of populations in twDEMC
	})

setMethodS3("getNChains","twDEMC", function( 
		### Extracts the number of chains
		res	##<< object of class twDEMC
		,... 
	){
		# calcNGen.twDEMC
		##seealso<<   
		## \code{\link{subChains.twDEMC}}
		## ,\code{\link{twDEMCInt}}
		ncol(res$rLogLik)
		### integer, number of populations in twDEMC
	})



setMethodS3("replacePops","twDEMC", function( 
		### replaces several populations by toher ones
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
		nChains = ncol(x$rLogLik)
		nChainsPop = nChains %/% nPop
		if( length(iPops) != nPopNew) stop("replacePops.twDEMC: replacing population has number of populations that differs from length of iPop.")
		nChainsPopNew = ncol(xNew$rLogLik) %/% nPopNew
		if( nChainsPop != nChainsPopNew) stop("replacePops.twDEMC: both populations must have the same number of chains per populations.")
		xN <- if( x$thin != xNew$thin ) thin(xNew,newThin=x$thin) else xNew
		if( !identical( dim(x)[1:2], dim(xN)[1:2]) ) stop("replacePops.twDEMC: both populations must have the same dimensions of variables and cases.")

		res <- x
		iChains <- matrix(1:nChains, nrow=nChainsPop)[,iPops]
		res$temp[,iPops] <- xN$temp
		res$parms[,,iChains] <- xN$parms
		res$rLogLik[,iChains] <- xN$rLogLik 
		res$pAccept[,iChains] <- xN$pAccept
		if( 0 < length(x$resFLogLikX))
			res$resFLogLikX[,iChains] <- xN$resFLogLikX
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
			window(	mcmc( t(adrop(x$parms[,,i,drop=FALSE],3)), thin=origThin ), start=start, thin=thin )
		})	
	mcmc.list( tmp.l )
	### a \code{\link[coda]{mcmc.list}}
}
#)
#mtrace(as.mcmc.list.twDEMC)


setMethodS3("subset","twDEMC", function( 
	### Condenses an twDEMC List to the cases boKeep.
	x, 
	boKeep,
	...
){
	# subset.twDEMC 
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	if( is.logical(boKeep) )
		boKeep <- rep(boKeep,length.out=nrow(x$rLogLik) )
	if( 0 < length(x$resFLogLikX)){
		#if dropped last row, set resFLogLikX to -Inf
		if( is.numeric(boKeep) )
			if( !(nrow(x$rLogLik)%in%boKeep) )
				x$resFLogLikX[] <- -Inf
		if( is.logical(boKeep) )
			if( !boKeep[nrow(x$rLogLik)] )
				x$resFLogLikX[] <- -Inf
	}
	x$parms <- x$parms[,boKeep,, drop=FALSE] 
	x$rLogLik <- x$rLogLik[boKeep,, drop=FALSE] 
	x$pAccept <- x$pAccept[boKeep,, drop=FALSE]
	if( !is.null(x$temp))
		if( is.null(ncol(x$temp)) )
			x$temp <- x$temp[boKeep, drop=FALSE]
		else
			x$temp <- x$temp[boKeep,, drop=FALSE]
	#better keep this attribute #attr(res,"batchCall") <- NULL
	x
	### list of class twDEMC with subset of cases
})
#mtrace(subset.twDEMC)

setMethodS3("thin","twDEMC", function( 
	### Reduces the rows of an twDEMC object (list returned by \code{\link{twDEMCInt}}) to correspond to a thinning of \code{newThin}.
	x, ##<< the twDEMC list to thin 
	newThin=x$thin, ##<< the target thinning factor, must be positive multiple of vMcpl$thin 
	start=1, ##<< the start time of the chain
	end=NULL, ##<< the maximum end time of the chains
	...
	, doKeepBatchCall=FALSE	##<< wheter to retain the batch call attribute of x

){
	# thin.twDEMC
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
		
	# with the thinned list having mZ rows, this corresponds to (mZ-1)*thin Metropolis steps + 1 row for the initial state
	if( (newThin < x$thin) | (newThin %% x$thin) )
		stop(paste("increased thin must be a multiple of former thin",x$thin))
	#thin own past: keep first line and every occurenc of multiple thin
	thinFac <- newThin %/% x$thin
	iGen1 <- (2:ncol(x$parms))*x$thin	#will keep the first row
	nGen1 <- length(iGen1)
	if( start > 1)
		keep.s <- iGen1 >= start 
	else keep.s <- rep(TRUE, nGen1)
	if( !is.null(end) )
		keep.s[iGen1 > end] <- FALSE
	#which( 1:5 %% 2 == 0, arr.ind = TRUE)+1
	keep.i <- 1+which( keep.s & (1:(ncol(x$parms)-1) %% thinFac == 0), arr.ind = TRUE )
	if( start == 1) keep.i <- c(1, keep.i)	#keep starting value
	res <- subset.twDEMC( x, keep.i )
	if( !keep.s[nGen1] & 0<length(res$resFLogLikX) ) res$resFLogLikX[] <- -Inf	#do not delete but keep names, but make sure that is accepted next time
	res$thin <- newThin
	if(!doKeepBatchCall) attr(res,"batchCall") <- NULL
	res
})
#mtrace(thin.twDEMC)



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
	res$rLogLik <- abind( lapply(.pops, function(psr){psr$rLogLik}) )
	if( 0 < length(x$resFLogLikX))
	res$resFLogLikX <- abind( lapply(.pops, function(psr){psr$resFLogLikX}) )
	res$pAccept <- abind( lapply(.pops, function(psr){psr$pAccept}) )
	res$thin <- x$thin
	res$temp <- abind( lapply(.pops, function(psr){psr$temp}), along=2 )
	for( i in c("parms","rLogLik","pAccept","temp") )
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
	x,
	...
){
	#stackChains.twDEMC
	##seealso<<   
	## \code{\link{subChains.twDEMC}}
	cbind( rLogLik = abind( lapply( 1:dim(x$rLogLik)[2], function(i){ x$rLogLik[,i] }), along=1 )
		,parms = stackChains.array(x$parms)
	)
	### Matrix with first column the log-Likelihold rLogLik and the remaining columns the variables.
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



twDEMCPopMeans <- function(
	### Calculating population means across chains, and smooth time series.
	x				##<< a matrix with columns chains
	,nPops			##<< number of populations 
	,kSmooth=NULL	##<< weights to the filter function, or just number of points
){
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
	if( boMat) xPop else as.vector(xPop)		
}

twDEMCPopApply <- function(
	### Applying a function across all chains of one population for each case.
	x				##<< a matrix with columns chains or array with last dimension chain
	,nPops			##<< number of populations
	,FUN			##<< function to apply to population submatrix
	,...			##<< further arguemtns to FUN
){
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
	### array with last dimenstion correponding to population
}
#twDEMCPopApply( res$rLogLik, nPop, apply, 1, mean )	#rLoglik for each case by population
#twDEMCPopApply( res$rLogLik, nPop, as.vector )			#stack rLogLik
#twDEMCPopApply( res$parms, nPop, function(x){ abind(twListArrDim(x),along=2) })	#stack param columns by population


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

getDiffLogLik.twDEMCProps <- function(
	### Extract the Differences in Log-Likelihood between accepted states and proposals.
	Y					##<< matrix of proposals with row "accepted" and first step (column) initial state, rows: results components of fLogLik and third dimension chains.
	,resCols			##<< the rows of Y with result components, either names or positions
	,nLastSteps = 128	##<< number of last steps of Y for which to extract diffs
	,temp=1				##<< numeric matrix (resComp x pop): the temperature applied to the difference of LogLiks
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
	acceptedPos <- getAcceptedPos.twDEMCProps(YL,...) 
	diffLogLik <- abind(lapply( 1:dim(YL)[3], function(j){ 	adrop((YL[resCols,-1,j,drop=FALSE] - YL[resCols,acceptedPos[-1,j],j,drop=FALSE])/temp[,(j-1)%/%nPopsChain+1],3)}),rev.along=0)
	diffLogLik
	### numeric array ( component x nLastSteps x chain ) of Lp-La, the L
}

replaceNonFiniteDiffLogliks <- function(
	### For each component replace NAs by sample of others and non-Finite values by minimum of others  
	diffLogLik			##<< array( streams x steps) Lp-La see \code{\link{getDiffLogLik.twDEMCProps}}
	,doConstrainNeg=FALSE	##<< if given, likelihood of accepted jumps (positive) is constrained to 0
){
	d <- t(apply( diffLogLik,1,function(ds){
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
	dimnames(d) <- dimnames(diffLogLik)
	if( doConstrainNeg )
		d[d>0] <- 0
	d
}

sampleAcceptedFixedTempDiffLogliks <- function(
	### Calculate diffLogLik for fixed temperature components and return subset of accepted ones 
	diffLogLik			##<< array( streams x steps) Lp-La see \code{\link{getDiffLogLik.twDEMCProps}}
	,TFix=numeric(0)
){
	dTFix <- colSums( diffLogLik[names(TFix),,drop=FALSE]/TFix )
	acceptedFixed <- dTFix>log(runif(ncol(diffLogLik)))
	#pa <- sum(acceptedFixed)/nj
	diffLogLik[,acceptedFixed,drop=FALSE]
}

