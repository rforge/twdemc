fDenMultiMode <- function(x){
	y <- +2.5*cos(4*x+pi/2) + 5*cos(2*x) + 9*cos(x)
	z <- -20*(y +  (x/pi)^4) -200
}

.tmp.test <- function(){
	x1 <- seq(-2*pi,+2*pi,length.out=300)
	#x1 <- seq(-10,+10,length.out=300)
	plot( fDenMultiMode(x1) ~ x1,  type="l" )
	x2 <- seq(-4.8,-1,length.out=300)
	plot( fDenMultiMode(x2) ~ x2,  type="l" )
}

.tmp.fit <- function(){
	thetaPrior = c(x=0)
	covarTheta = c(x=1/pi)
	dInfos = list(d1=list(fLogDen=fDenMultiMode))
	nObs=1
	#
	gelmanCrit = 1.4
	dLogDenCrit = 2		##<< difference in median logDensity above which pops are regarded significantly different
	set.seed(0816)
	res <- res4 <- twDEMCSA( thetaPrior=thetaPrior, covarTheta=covarTheta, nPop=2, dInfos=dInfos, nObs=nObs
	, nGen=32*5, ctrlBatch=list(nGen0=32,nGenBatch=32)
	, ctrlConvergence = list(gelmanCrit=gelmanCrit	, critSpecVarRatio=80	)
	)
mcs <- concatPops(stackChainsPop(stackPopsInSpace(res)))
plot( as.mcmc.list(mcs))
gelman.diag(as.mcmc.list(mcs))

	resTodo <- list(res0=res)
	while( length(resTodo) ){
		res1 <- resTodo[[ length(resTodo) ]]
		res <- twDEMCSA( res1, nGen=32 )
		resTodo[ length(resTodo) ] <- NULL
		#	
		mcs <- concatPops(stackPopsInSpace(stackChainsPop(res)))
		gelmanDiagRes <- try( {tmp<-gelman.diag(mcs); if(length(tmp$mpsrf)) tmp$mpsrf else tmp$psrf[1]} )	# cholesky decomposition may throw errors
		if( gelamnDiagRes > gelmanCrit){
			# check if temperated LogDen differs
			logDenT <- calcTemperatedLogDenChains(mcs)
			tmp <- apply(logDenT,c(1,3),sum)
			tmp2 <- tmp[ apply(tmp,1,min) > quantile(tmp,0.2), ]
			matplot(tmp2, type="l")
			ldm <- apply(tmp,2,median)
			oldm <- order(ldm, decreasing=TRUE)	# highest first
			spacesPopUnique <- unique(getSpacesPop(res))
			if( diff(ldm) > dLogDenCrit ){
				# store other one in Todo-List
				resTodo[[as.character(length(resTodo)+1)]] <- 
					.splitPopsSpace(res, iSpace=spacesPopUnique[ oldm[2] ]) 
			}
			# store the best one at the end of the Todo list (so continue with that)
			resTodo[[as.character(length(resTodo)+1)]] <- 
				.splitPopsSpace(res, iSpace=spacesPopUnique[ oldm[1] ]) 
		} else{
			# check logDenDrift not valid when chains or pops are stacked and time is lost
			
			# Todo: .stackChainsWithinPop
			
			tld1 <- rowSums(calcTemperatedLogDen(subsetTail(res1,0.1)))
			tld2 <- rowSums(calcTemperatedLogDen(subsetTail(res,0.1)))
			boxplot( tld1,tld2 )
			if( isLogDenDrift(logDenT[,1,], res$dInfos) )
		
			# check shift in logDen
			isLogDenDrift
			if( )
			apply(tmp,2,quantile, probs=c(0.025,0.975))	# despite several modes, same spread of logDen
			
			.spacesPop <- getSpacesPop(res)
			iSpaces=.spacesPop[c(1)]
			#iSpaces=.spacesPop[c(2)]
			resSpaces <- subPops(res, iSpace=c(iSpaces))
			# split the population into two independent pops 
			fKeepOdd <- function(pop){ (seq(1,nrow(pop$parms)) %% 2) == 1 } 
			fKeepEven <- function(pop){ (seq(1,nrow(pop$parms)) %% 2) == 0 } 
			tmp1 <- subsetF( resSpaces, fKeepOdd )
			tmp2 <- subsetF( resSpaces, fKeepEven )
			tmp2$pops <- lapply(tmp2$pops, function(pop){ pop$spaceInd <- 5; pop})
			resSpaces$pops <- c( tmp1$pops, tmp2$pops )
			
		} # if gelmanDiag
	} # while Pops
}

.splitPopsSpace <- function( 
	### split all the the population into two independent pops
	mc			##<< the twDEMCPops object for which to split a subspace
	,iSpace=spacesPopUnique[1]	##<< the subSpace which should be splitted
){
	spacesPopUnique <- unique(getSpacesPop(mc))
	mcSpace <- subPops(mc, iSpace=iSpace)
	fKeepOdd <- function(pop){ (seq(1,nrow(pop$parms)) %% 2) == 1 } 
	fKeepEven <- function(pop){ (seq(1,nrow(pop$parms)) %% 2) == 0 } 
	tmp1 <- subsetF( mcSpace, fKeepOdd )
	tmp2 <- subsetF( mcSpace, fKeepEven )
	tmp2$pops <- lapply(tmp2$pops, function(pop){ pop$spaceInd <- max(spacesPopUnique)+1; pop})
	mcSplit <- mcSpace
	mcSplit$pops <- c( tmp1$pops, tmp2$pops )
	mcSplit
}

