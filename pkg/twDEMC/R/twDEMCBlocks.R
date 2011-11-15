twDEMCBlockInt <- function(
	### Differential Evolution Markov Chain with blocked parameter update
	pops = list( list( ## list of population infos for each population, each info is list with components
		##describe<<
		Zinit		##<< list of matrices (nState x nParm  x nChain) initial states for each population see details and \code{\link{initZtwDEMCNormal}}.
		, nGen = 10	##<< number of generations, i.e. steps to proceed
		, T0=1		##<< initial temperature
		, Tend=1	##<< end temperature
		, Tprop=NULL	##<< temperature proportions of result components
		, X=NULL		##<< numeric matrix (nParm x nChain) initial state
		, logDenCompX=NULL 		##<< numeric matrix (nComp x nChain): logDen components of initial state X, see details
	)) ##end<<
	, blocks = list( list( ##<< list of parameter blocks, each block is list with entries
			##describe<<
			, compPos=1:nrow(pops[[1]]$Zinit)	##<< index or names of the parameter components to be updated
			, fLogDen=NULL
				###	\code{function(theta, ...)} calculates a vector of logDensities 
				### corresponding to different data streams of parameter vector theta 
				### \cr or \code{function(theta, logDenCompX, metropolisStepTemp, ...)} 
				### to handle first steps in Multi-step Metropolis decision internally. 
				### See details.  
			, fUpdateBlock=.doDEMCStepsBlock	##<< function to update the parameters
			, argsFLogDen=list()	##<< further arguments passed to fLogDen
			, argsFUpdate=list()	##<< further arguments passed to fUpdate
			, isNeedSteps=TRUE		##<< TRUE, if jumps for Metropolis proposals should be generated for this block
			, intResComp=vector("integer",0) 
				### integer or character vector: indices or names of results components of fLogDen 
				### that are used for internal Metropolis decisions
			, fLogDenScale=1 
				### scalar multiplied to the result of fLogDen 
				### allows using functions of negative LogDensity (-1) or Gaussian misfit function (-1/2) instead of logDensity
			, TFix = vector("numeric",0) ##<< named numeric vector with Temperature for result components that have fixed temperature
		)) ##end<<
	, nGen = 10			##<< number of generations, default, if not given in pops[[i]] description
	, controlTwDEMC = list()	##<< DEMCzsControl parameters influencing the update and the convergens of the chains (see details)	
	, debugSequential=FALSE 		##<< if TRUE apply is used instead of sfApply, for easier debugging
	, remoteDumpfileBasename=NULL	##<< the basename of a dumpfile that is created on error on remote process 
	, fDiscrProp=NULL				##<< function applied to proposal, e.g. to round proposals to to discrete possible values function(theta,...)
	, argsFDiscrProp=list()			##<< further arguments to fDiscrProp
	, doRecordProposals=FALSE		##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y
	, upperParBounds = vector("list",nPops)
	### list of named numeric vectors: giving upper parameter bounds for each population 
	### for exploring subspaces of the limiting distribution, see details
	, lowerParBounds = vector("list",nPops)  ##<< similar to upperParBounds
){
	ctrl = list( 
		thin = 4 	   	##<< thinning interval	 
		#moved to block,TFix = vector("numeric",0)		##<< named numeric vector: result components for which temperature shoudl not change		
		,F = 2.38 		##<< related to multiplicative error (F2=F/sqrt(2*Npar), see eps.mult 
		,pSnooker= 0.1,	##<< probability of a snooker update (others parallel updates)
		pGamma1 = 0.1,	##<< probability of jumping to state of another chain (different modes)
		epsMult =0.2,	##<< >0 gives d-dimensional variation around gamma. It adds scaled uncorrelated noise to the proposal. Its advantage over eps.add is that its effect scales with the differences of vectors in the population whereas eps.add does not. if the variance of a dimensions is close to 0, eps.mult gives smaller changes. \cr A uniformly distributed error, i.e. F2*runif(1+-epsMult*prop) multiplied to difference vector from parallel update 
		epsAdd = 0,	   	##<< >0 is needed to ensure that all positions in the space can be reached. For targets without gaps, it can set small or even to 0. \cr sd of normally distributed error added to proposal by parallel or snooker update. 
		pAcceptWindowWidth = 256, ##<< number of generations back over which the acceptance rate is calculated
		probUpDir=0.5 	##<< probability of direction between two states of increasing Density (increase during burin may accelerate convergence)
		,initialAcceptanceRate=0.25	##<< numeric vector (nChains) initially assumed acceptance rate. Used to calculate the number of generations backwards to sample from
		,DRgamma=0		##<< factor for reducing step length [0..1) in delayed rejection step, 0 means no DR step
		,minPCompAcceptTempDecr=0.15  ##<< if acceptance rate drops below minPCompAcceptTempDecr+0.02 times this level, employ delayed rejection (DR)
		,pIndStep = 1.5 ##<< independent state about after on average about those number of 1.5 accepted steps
		,nPastGen = 10  ##<< factor for determining the number of recent past states to sample during burnin. It is multiplied by the number of parameters. Past generations are calculated by deviding by the number of chains per population 
	)  ##end<< 
	ctrl[names(controlTwDEMC)] <- controlTwDEMC

	#--  dimensions of Zinit
	pops <- lapply( pops, .checkPop, nGenDefault=nGen) # fill in default values for missing entries 
	ZinitPops <- lapply(pops,"[[","Zinit")
	parNames <- colnames(ZinitPops[[1]])
	nParm <- ncol(ZinitPops[[1]])
	nPop <- length(pops)
	nChainPop <- dim(ZinitPops[[1]])[3]
	nChain <- nChainPop*nPop 
	iPops <- 1:nPop
	popChain <- rep(1:nPop,each=nChainPop)	# population for chain at given index
	chainsPop <- lapply( iPops, function(iPop){ (iPop-1)*nChainPop+(1:nChainPop)}) # chains for given population
	M0Pops <- sapply( ZinitPops, nrow )
	nGenPops <- sapply( pops, "[[", "nGen")
	nThinnedGenPops = sapply( nGenPops %/% ctrl$thin, max, 1 )	#number of thinning intervals 
	nGenPops = nThinnedGenPops * ctrl$thin  #number of generations in this batch (do not need to move behind last thinning interval)
	Nz <- M0Pops +(nThinnedGenPops)			   #number of rows in Z after batch
	XPops <- lapply( pops, "[[", "X" )
	XChains <- do.call(cbind,XPops)
	logDenCompXPops <- lapply( pops, "[[", "logDenCompX" ) # see initial states for initializing missing entries
	
	#-- dimensions of blocks
	blocks <- lapply( blocks, .checkBlock, parNames=parNames)
	nBlock <- length(blocks)
	iBlocks <- 1:length(blocks)
	#block <- blocks[[1]]
	compPosBlock <- lapply( blocks, "[[", "compPos" ) # already transformed to positions in .checkBlock
	if( (length(compPosBlock) != nParm) || !all.equal( sort(do.call(c,compPosBlock)), 1:nParm, check.attributes = FALSE) ) 
		stop("union of all blocks must equal parameters (columns of Zinit)")
	
	#-- initialize further parameters to parallel and snooker steps
	ctrl$Npar12  =(nParm - 1)/2  # factor for Metropolis ratio DE Snooker update
	ctrl$F2 = ctrl$F/sqrt(2*nParm)
	ctrl$F1 = 1
	ctrl$gInd <- ctrl$nPastGen*nParm/nChainPop	#number of independent generations to sample from, similarly to M0 (usually ctrl$nPastGen=10)
	# terBraak report that may not sample the distribution, if not using the full past
	# but together with decreasing temperature acceptance rate drops very low
	# hence constrain to the past during burnin
	

	#-- evaluate logDen of initial states and get names of the result components
	##details<< \describe{ \item{Initial state: \code{X}}{
	## If initial state X is not specified, the last column (generation) of Z is utilized.
	## If in addition to X, logDenX is specified as a numeric vector, the fLogDen will not be avaluated for the initial state of the chains.
	## All the results of fLogDen for the initial state must be finite.
	## }}
	# find those chains where no logDenCompX is given, or not all components are finite
	popsMissingResX <- ( 0 == sapply( logDenCompXPops, length) )
	iMissingPops <- which(popsMissingResX)
	iNonMissingPops <- which(!popsMissingResX)
	missingPopInfo <- do.call( rbind, lapply( iMissingPops, function(iPop){
				cbind(iPop=iPop,iChainInPop=1:nChainPop)
			} ))
	missingPopInfo <- rbind( missingPopInfo, do.call( rbind, lapply( iNonMissingPops, function(iPop){
				chainIsAllFinite <- apply( logDenCompXPops[[iPop]], 2, all, is.finite )
				cbind(iPop=iPop,iChainInPop=which(chainIsAllFinite))
			})))
	if( 0 < nrow(missingPopInfo) ){
		# evaluate logDen for those
		missingPopInfo <- cbind(missingPopInfo, iChain=(missingPopInfo[,"iPop"]-1)*nChainPop + missingPopInfo[,"iChainInPop"])
		XChainsMissing <- XChains[,missingPopInfo[,"iChain"] ]	# initial states of the missings
		logDenCompBlock <- vector("list", nBlock )		
		for( iBlock in iBlocks ){
			block = blocks[[iBlock]]
			#tmp <- block$fLogDen
			#tmp2 <- do.call( tmp, c(list(XChainsMissing[,1]),block$argsFLogDen))
			#mtrace(twCalcLogDenPar)
			.resLogDenPar <- twCalcLogDenPar(
				fLogDen=block$fLogDen
				, xProp=t(XChainsMissing[,,drop=FALSE])
				, logDenCompX=numeric(0)
				, argsFLogDen=block$argsFLogDen
				, fLogDenScale=block$fLogDenScale	
				, intResCompNames=block$intResCompNames	
				, debugSequential=debugSequential
				, remoteDumpfileBasename=remoteDumpfileBasename 
			)
			if( any(!is.finite(.resLogDenPar$logDen)) )
				stop(paste("twDEMCBlockInt: non-finite logDensity of starting value ",sep=""))
			blocks[[iBlock]]$resCompNames <- .getResFLogDenNames( t(.resLogDenPar$logDenComp) ) 
			logDenCompBlock[[iBlock]] <- .resLogDenPar$logDenComp
		} #iBlock
		logDenCompXMiss <- do.call( cbind, logDenCompBlock )
		for( jPop in seq_along(iMissingPops) ){
			iPop <- iMissingPops[jPop]
			iInfo <- which(  missingPopInfo[,"iPop"] == iPop )
			logDenCompXPops[[ iMissingPops[jPop] ]] <- 
				t(logDenCompXMiss[iInfo,])
		}
		for( jPop in seq_along(iNonMissingPops) ){
			iPop <- iNonMissingPops[jPop]
			iInfo <- which(  missingPopInfo[,"iPop"] == iPop )
			logDenCompXPops[[iPop]][,missingPopInfo$iChainInPop[iInfo] ] <-
				t(logDenCompXMiss[iInfo,])
		}
	}else{ # any isMissing
		# no Missings, evaluate fLogDen to obtain result component names
		for( iBlock in iBlocks ){
			block = blocks[[iBlock]]
			.resFLogDen <- do.call( block$resFLogDen, c(list(XChains[,1]),block$argsFLogDen))
			blocks[[iBlock]]$resCompNames <- .getResFLogDenNames( .resFLogDen ) 
		} #iBlock
	}
	logDenCompXChains <- do.call( cbind, logDenCompXPops )
	
	#-- get indices of internal components
	for( iBlock in iBlocks ){
		block = blocks[[iBlock]]
		intResCompBlockI <- block$intResComp
		if( is.character(intResCompBlockI) ){
			blocks[[iBlock]]$intResComp <- match(intResCompBlockI, block$resCompNames )
		}
	} #iBlock

	#-- make unique resCompNames and calculate positions in resCompNames of results of different blocks
	resCompNamesOrig <- lapply( blocks, "[[", "resCompNames" )
	resCompNamesUnique <- lapply( iBlocks, function(iBlock){ 
			.getRescompNameBlock(resCompNamesOrig[[iBlock]], iBlock)
		})
	resCompNamesFlat <- do.call(c, resCompNamesUnique)	# concatenate names
	nResComp <- length(resCompNamesFlat)
	resCompNamesPos <- { # list: for each block: position of resCompNames in resCompNamesFlat
		tmp.length <- sapply(resCompNamesUnique, length)
		tmp.start <- cumsum(tmp.length)-tmp.length[1]
		tmp.i <- lapply(iBlocks,function(iBlock){ tmp.start[iBlock]+1:tmp.length[iBlock] })
	}
		
	
	#-- calculate temperature steps: exponentially decreasing from T0 to Tend
	# temperature for T0 and then for each (unthinned) generation
	temp <- lapply( iPops, function(iPop){
			c( pops[[iPop]]$T0, pmax(1,calcDEMCTemp( pops[[iPop]]$T0, pops[[iPop]]$Tend, nGenPops[iPop] )) )			
		})
	
	#-- setup acceptance rate recording
	# number of requried independent generations to choose from (10*d independent states: TerBraak06 TerBraak08) devided by number of chains
	if( length(ctrl$initialAcceptanceRate) == 1)
		ctrl$initialAcceptanceRate <- matrix(ctrl$initialAcceptanceRate, nrow=nBlock, ncol=nChain, dimnames=list(blocks=NULL, chains=NULL) )
	if( !is.matrix(ctrl$initialAcceptanceRate) || dim(ctrl$initialAcceptanceRate) != c(nBlock,nChain))
		stop("ctrl$initialAcceptanceRate must be matrix of dim nBlock x nChain")
	arPops <- apply(popMeansTwDEMC(ctrl$initialAcceptanceRate, nPop),1,pmax, 1/ctrl$thin) # current acceptance rate of population
	#nGenBack <- pmin(mZ,ceiling(ctrl$rIndToAcc * ctrl$gInd * pmax(1,ctrl$pIndStep/(ar * ctrl$thin))))	#number of genrations to select states from for each population, ctrl$gInd is multiplied by number of rows for one step depending on acceptance rate and thinning but at least one
	arMinPops <- apply(arPops,2,min)	# minimum of acceptance rate across blocks for each population
	nGenBackPops <- pmin(M0Pops,ceiling(ctrl$gInd * pmax(1,ctrl$pIndStep/(arMinPops * ctrl$thin))))	#number of genrations to select states from for each population, ctrl$gInd is multiplied by number of rows for one step depending on acceptance rate and thinning but at least one  
	##details<< \describe{ \item{Acceptance rate}{
	## The acceptance rate is tracked for each chain across ctrl$pAcceptWindowWidth generations. \cr
	## If acceptance rates drops to low values, this might be because of bad localization
	## ,i.e the requirement to use states from more distant past.
	## In order to improve localization, less parameters or more chains per population are required.
	## }}
	# acceptWindow holds the number of accepted steps for each thinning interval for each chain
	# if its boundaries are filled up, the last nGenBack states are copied to the first interval
	aThinWinW <- ctrl$pAcceptWindowWidth %/% ctrl$thin	# number of thinning intervals over which to calculate average acceptance rate
	ctrl$pAcceptWindowWidth <- aThinWinW * ctrl$thin
	# acceptWindow records number of accepted steps in thinnging interval
	#acceptWindow <- matrix( rep(ctrl$thin*ar, each=2*aThinWinW*nChainsPop), nrow=2*aThinWinW, ncol=d$chains )	#record number of accepted steps in thinning interval
	acceptWindow <- array( rep(ctrl$thin*ctrl$initialAcceptanceRate, each=2*aThinWinW)
	, dim=c(2*aThinWinW, dim(ctrl$initialAcceptanceRate) ) 
	, dimnames=c(list(steps=NULL),dimnames(ctrl$initialAcceptanceRate)) )	#record number of accepted steps in thinning interval
	
	#-- preallocate output of parameters, calculated LogDen, indices of accepted steps, and next proposal
	# need lists for populations because Nz (number of samples including intial population ) and nGen (samlples without initial) and nThinnedGen is differ across populations in general	
	ZPops <- lapply( iPops, function(iPop){
			dNames <- if( !is.null(dimnames(ZinitPops[[iPop]])) ) dimnames(ZinitPops[[iPop]]) else list( samples=NULL, parms=paste("p",1:nParm,sep="_"), chains=NULL )
			Z <- array(NA_real_, dim=c(Nz[iPop],nParm,nChainPop), dimnames=dNames)
			Z[1:M0Pops[iPop],,] <- ZinitPops[[iPop]]
			Z
		})
	# several resLogDen per chain (blocks are concatenated)
	resLogDen = lapply(iPops,function(iPop){ res <- array( NA_real_
				, dim=c(1+nThinnedGenPops[iPop], nResComp, nChainPop)
				, dimnames=list(stepsThin=NULL, resComp=resCompNamesFlat, chains=NULL))
			res[1,,] <- logDenCompXPops[[iPop]]
			res
		}) 
	# one rLogDen per chain and per block
	logDen <- lapply(iPops, function(iPop){	res <- array(NA_real_
				, dim=c(1+nThinnedGenPops[iPop], nBlock, nChainPop) 
				, dimnames=list(stepsThin=NULL, block=NULL, chains=NULL) )
			for( iChainPop in 1:nChainPop)
				for( iBlock in 1:nBlock )
					res[ 1,iBlock,iChainPop] <- sum(logDenCompXPops[[iPop]][resCompNamesPos[[iBlock]],iChainPop ])
			res
		})
	pAccept <- logDen
	for( iPop in iPops ){
		pAccept[[iPop]][1,,] <- ctrl$initialAcceptanceRate[,chainsPop[[iPop]] ]
	}
	#pAccept[1,,] <- NA_real_ 	# TODO: initialize acceptance rate
	# record of proposals and fLogDen results, rows c(accepted, parNames, resCompNames, rLogDen)
	# if doRecordProposals is FALSE record only the thinning intervals for the last 128 generations
	# +1 for the initial state
	nThinLastPops <- if(doRecordProposals) nThinnedGenPops else pmin(nThinnedGenPops, ceiling(128/ctrl$thin))
	nThinOmitRecordPops = nThinnedGenPops-nThinLastPops	#the Thinning intervals with no recording of outputs
	nGenOmitRecordPops = nThinOmitRecordPops*ctrl$thin
	nGenYPops <- nThinLastPops*ctrl$thin #+1?
	YPops <- lapply(iPops, function(iPop){ res <- array( NA_real_ 
				, dim=c( nGenYPops[iPop], nParm+nBlock+nResComp, nChainPop)
				, dimnames=list(steps=NULL, vars=c(rownames(XPops[[iPop]]),paste("accepted",iBlocks,sep=""),colnames(resLogDen[[iPop]])), chains=NULL ) )
			if( nGenYPops[iPop] == nGenPops[iPop] ) # record first state as accepted
				res[1,,] <- rbind( XPops[[iPop]], matrix(TRUE, nrow=nBlock, ncol=nChainPop), logDenCompXPops[[iPop]] )
			res
		})
	
	
	
	#-- return
	##value<< a list of populations, each entry is a list
	res <- lapply( iPops, function(iPop){ list( 
			##describe<<
			Z = ZPops[[iPop]]			##<< numeric array (steps x parms x chains): collected states, including the initial states
			,temp = temp[[iPop]]		##<< numeric vector: global temperature, i.e. cost reduction factor
			,pAccept= pAccept[[iPop]]	##<< acceptance rate of chains
			,resLogDen = resLogDen[[iPop]]	##<< numeric array (steps x resComps x chains): results components of fLogDen of blocks  
			,logDen = logDen[[iPop]]	##<< numberic array (steps x block x chains): results summed over blocks
			,Y = YPops[[iPop]]
		)}) ##end<<
}
attr(twDEMCBlockInt,"ex") <- function(){
	data(twLinreg1); attach( twLinreg1 ) 
	
	.nPops=2
	Zinit <- ZinitPops <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=4*.nPops, nPops=.nPops)
	pops <- list(
		pop1 <- list(
			Zinit = ZinitPops[1:3,,1:4,drop=FALSE]	# the first population with less initial conditions
			,nGen=10
		),
		pop2 <- list(
			Zinit = ZinitPops[,,5:8,drop=FALSE]	# the first population with less initial conditions
			,nGen=15
			,T0=10
		)
	)
	#tmp <- .checkPop(pops[[1]])
	
	# both blocks compare the same model against the same observations
	blockDefault <- list(
		fLogDen=logDenGaussian
		,resCompNames=c("obs","parms")
		,TFix=c(parms=1)
		,argsFLogDen <- list(
			fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
			obs=obs,				### vector of data to compare with
			invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
			xval=xval
		)
	)
	# first block updates only parameter a, second block updates only parameter b
	# contrary to this example, correlated parameters should not be distributed across blocks
	blocks <- list(
		blockA = within(blockDefault,{
				compPos="a"
				argsFLogDen <- c( argsFLogDen, list(
						thetaPrior= thetaTrue["a"]	### the prior estimate of the parameters
						,invCovarTheta = invCovarTheta[1,1,drop=FALSE]	### the inverse of the Covariance of the prior parameter estimates
						,blockIndices=1
					)) 
			})
		,blockB = within(blockDefault,{
			compPos="b"
			argsFLogDen <- c( argsFLogDen, list(
					thetaPrior= thetaTrue["b"]	### the prior estimate of the parameters
					,invCovarTheta = invCovarTheta[2,2,drop=FALSE]	### the inverse of the Covariance of the prior parameter estimates
					,blockIndices=2
			)) 
		})
	)
	#tmp <- .checkBlock(blocks[[2]],parNames=c("a","b"))
	block <- blocks[[2]] 
	block <- blocks[[1]] 
	#mtrace(logDenGaussian)
	(tmp <- do.call( logDenGaussian, c(list(pops[[1]]$Zinit[1,,1]), block$argsFLogDen) ))
	
	
	#mtrace(twDEMCBlockInt)
	res <- twDEMCBlockInt( pops=pops, blocks=blocks, nGen=60)
	str(res[[2]])
	#plot(res[[2]]$temp)
	
	# load the restart file and continue
	load( file=restartFilename )	# variable resRestart.twDEMC
	getNGen(resRestart.twDEMC)		# 48, last thinnging interval (each thin=4) before 50
	res2 <- twDEMCBatch(resRestart.twDEMC)	# continue without needing to respecify parameters
	getNGen(res2)					# 60 as nGen
	res3 <- twDEMCBatch(res2, nGen=100)	    # continue even further without needing to respecify parameters
	getNGen(res3)					# 100 as nGen
	
	detach()
	
}

.getRescompNameBlock <- function(
	### appending iBlock to the resCompNames to make the names unique among blocks
	resCompNames
	,iBlock
){
	paste(resCompNames,iBlock,sep="_")
}

.checkPop <- function(
	### filling default values in pop description
	pop
	,nGenDefault=10
){
	##details<<
	## If several entries are null, they are set to default values
	if( 0==length(pop$Zinit) || !is.numeric(pop$Zinit) || length(dim(pop$Zinit)) != 3 )
		stop(".checkPop: entry Zinit must be a numeric 3D array")
	if( 0==length(pop$nGen) ) pop$nGen <- nGenDefault
	if( 0==length(pop$T0) ) pop$T0 <- 1
	if( 0==length(pop$Tend) ) pop$Tend <- 1
	if( 0==length(pop$X) ){
		pop$X <- pop$Zinit[nrow(pop$Zinit),,]
	}else{
		if( !all.equal( rownames(pop$X), colnames(pop$Zinit), check.attributes = FALSE) )
			stop(".checkPop: rownames of X and colnames of Zinit must match")
	}
	### pop with entries nGen, T0, Tend and X defined 
	pop	
}

.checkBlock <- function(
	### filling default values in block description
	block
	,parNames	##<< variable names
){
	if( 0 == length(block$compPos) ) stop(".checkBlock: entry compPos must be specified, either as index or as variable names")
	if( 0==length(block$resCompNames) || !is.character(block$resCompNames) ) 
		stop(".checkBlock: must specify names of result components in entry resCompNames")
	cpOrig <- block$compPos
	if( is.character(cpOrig) ) block$compPos <- match(cpOrig,parNames)
	if( any(is.na(block$compPos)) ){
		stop(".checkBlock: some of compPos specified in block are not in varNames of initial values.")
	}
	if( 0==length(block$isNeedSteps) ) block$isNeedSteps <- TRUE
	if( 0==length(block$intResCompNames) ) block$intResCompNames <- vector("character",0)
	if( 0==length(block$fLogDenScale) ) block$fLogDenScale <- 1
	if( 0==length(block$TFix) ) block$TFix <- vector("numeric",0)
	### block with entries isNeedSteps, intResCompNames, fLogDenScale, and TFix defined
	block
}

.getResFLogDenNames <- function(
	### Extracting/Creating names for result components of logDenComp
	logDenComp		##<< vector or array (vars in rows) of results of logDensity function 
	,blockName=NULL
	,blockInd=1
){
	##details<< 
	## Names of the result components of fLogDen are used to ditinguish internal components.
	## If the density function does not provide names, they are created.
	## If the density function has only one result component, the name is lost during parallel execution.
	## Hence, a single name is also re-created, to avoid errors on checking.
	if( 0 < length(blockName) ) blockName <- paste("block",blockInd)	
	if( is.array(logDenComp)){
		res <- rownames(logDenComp)
		return( if( !is.null(res) ) res else paste(blockName,1:nrow(logDenComp),sep="") )
	}else{
		res <- names(logDenComp)
		return( if( !is.null(res) && length(res) > 1) res else paste(blockName,1:length(logDenComp),sep="") )
	}
}


updateBlockDEMC <- function(
	### updating parameter block by differential evolution 
){
	
}

