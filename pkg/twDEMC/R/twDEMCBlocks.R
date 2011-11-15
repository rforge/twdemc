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
			, fUpdateBlock=.updateBlockTwDEMC	##<< function to update the parameters
			, argsFLogDen=list()	##<< further arguments passed to fLogDen
			, argsFUpdate=list()	##<< further arguments passed to fUpdate
			, useMetropolis=TRUE	##<< TRUE, if jumps for Metropolis proposals should be generated for this block
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
	, upperParBoundsPop = vector("list",length(pops))
	### list of named numeric vectors: giving upper parameter bounds for each population 
	### for exploring subspaces of the limiting distribution, see details
	, lowerParBoundsPop = vector("list",length(pops))  ##<< similar to upperParBounds
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
	
	#-- get indices of internal components
	intResCompBlock <- vector("list",nBlock)
	for( iBlock in iBlocks ){
		block = blocks[[iBlock]]
		intResCompBlockI <- block$intResComp
		if( is.character(intResCompBlockI) ){
			resI <- blocks[[iBlock]]$intResComp <- match(intResCompBlockI, block$resCompNames )
			if( any(is.na(resI)))
				stop(paste("not all names of intResCompNames (",paste(intResCompBlockI,collapse=","),") in return of fLogDen: ",paste(block$resCompNames,collapse=",")))
		}
		intResCompBlock[[iBlock]] <- blocks[[iBlock]]$intResComp
	} #iBlock
	
	#-- proportions of temperature numeric matrix (comp x populations) for temperature varying across data streams
	# XXTODO
	tmp.f <- function(){
		Tprop=matrix(1, nrow=nrow(logDenCompX), ncol=nPop, dimnames=dimnames(logDenCompX))
		if( length(ctrl$Tprop) > 1){
			if( !is.matrix(ctrl$Tprop) )
				ctrl$Tprop <- matrix( ctrl$Tprop, nrow=nrow(logDenCompX), ncol=nPop, dimnames=dimnames(logDenCompX) )
			Tprop <- ctrl$Tprop[ .getResFLogDenNames(logDenCompX[,1]), ,drop=FALSE]
			if( nrow(Tprop) != nrow(logDenCompX))
				stop("ctrl$Tprop must have a named entry for each component of fLogDen")
			Tprop[] = apply(Tprop,2,function(Tprop){Tprop / max(Tprop)}) #scale to maximum 1 per population
		}
	}
	
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
	# deprecated: temperature for T0 and then 
	# for each (unthinned) generation
	temp <- lapply( iPops, function(iPop){
			#c( pops[[iPop]]$T0, pmax(1,calcDEMCTemp( pops[[iPop]]$T0, pops[[iPop]]$Tend, nGenPops[iPop] )) )			
			pmax(1,calcDEMCTemp( pops[[iPop]]$T0, pops[[iPop]]$Tend, nGenPops[iPop] )) 			
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
	
	#-- arguments to fUpdate argument argsFUpdate that do not change within thinning interval
	argsFUpdateDefault = list(
		ctrl=ctrl
		,fCalcComponentTemp=calcComponentTemp 
		,fDiscrProp=fDiscrProp
		,argsFDiscrProp=argsFDiscrProp
	)
	# construct an list for each block that includes compPos, TFix and intResComp
	argsFUpdateBlocks <- lapply( blocks, function(block){ 
			c( argsFUpdateDefault, block )
		})
	argsUpdateBlocksTwDEMC <- list( 
		remoteFun=.updateBlocksTwDEMC,
		argsUpdateBlocksTwDEMC=list(
			argsFUpdateBlocks = argsFUpdateBlocks	# doDEMCStep
			,upperParBoundsPop=upperParBoundsPop
			,lowerParBoundsPop=lowerParBoundsPop
			,temp=temp
		)
	)
	argsUpdateBlocksTwDEMC$remoteDumpfileBasename<-remoteDumpfileBasename	#if null delete
	if( !debugSequential & sfParallel() )
		sfExport("argsUpdateBlocksTwDEMC")	#evaluated in remote process, only passed one time
	#else( assign("argsUpdateBlocksTwDEMC", argsUpdateBlocksTwDEMC, pos=1))	#export to current 
	tmp.remoteFunArgs <- if( !debugSequential & sfParallel() ) as.name("argsUpdateBlocksTwDEMC") else argsUpdateBlocksTwDEMC 	# eval.parent fails for argsUpdateBlocksTwDEMC within sequential function
	
	# arguments to demc that change between thinning intervals
	# substract logDen components of logDenCompX from logDenXExt
	chainState <- list(
		X = XChains
		,logDenCompX = logDenCompXChains
		#,logDenX = logDenX
	)
	
	for( iThin0 in (0:(min(nThinnedGenPops)-1)) ){
		mZPops <- M0Pops + iThin0
		# sample random vectors (note here parameters in rows)
		zxl <- lapply( iPops, function(iPop){
				.sampleStates(Z=ZPops[[iPop]],mZ=mZPops[iPop],nGenBack=nGenBackPops[[iPop]], nSample=ctrl$thin )				
			})	
		zx <- abind( zxl, along= 2 )
		# calculate proposed steps (differences not destinations) within next thinning interval
		genPropRes <- .generateXPropThin(zx,ctrl=ctrl, nChain=nChain)
		iGen = iThin0*ctrl$thin+(1:ctrl$thin)
		.tmp.f <- function(){
			tempThinStepsL <- lapply( iPops, function(iPop){
					matrix( temp[[iPop]][ iGen ], nrow=length(iGen), ncol=nChainPop )
				}) 
			tempThinSteps <- t(do.call( cbind, tempThinStepsL))	# steps need to be in last dimension (so chains in first)
			#tempThinSteps = lapply( temp, "[", t(temp[iGen,rep(1:nPop,each=nChainsPop),drop=FALSE])	#chains must be first dimension in order to acces by temp[i]
		}
		iPopsRecord <- which(iThin0 == nThinOmitRecordPops)
		for( iPop in iPopsRecord ){
			# record the initial state of the proposals record
			YPops[[iPop]][1,seq_along(parNames),] <- chainState$X[ ,chainsPop[[iPop]] ]
			YPops[[iPop]][1,nParm+iBlocks,] <- 1	#TRUE
			YPops[[iPop]][1,nParm+nBlock+seq_along(resCompNamesFlat),] <- chainState$logDenCompX[ ,chainsPop[[iPop]] ]
		}
		boRecordProposalsIThinPop = (iThin0 >= nThinOmitRecordPops)	# only record the last proposals after nThinOmitRecord
		
		iStep <- 1
		iChain <- 1
		iPop <- popChain[iChain]
		Xc <- chainState$X[,1]
		step <- genPropRes$xStep[,iChain,iStep]
		rExtra <- genPropRes$rExtra[iChain,iStep]
		for( iBlock in iBlocks){
			block <- blocks[[iBlock]]
			argsFUpdateBlock <- c( argsFUpdateBlocks[[iBlock]], list(
					step=step
					, rExtra=rExtra
					, temp=temp[[iPop]][iStep]
					, TProp=1 	# XXTODO
					, pAccept=0.25 # XXTODO
					, upperParBounds=upperParBoundsPop[[iPop]]
					, lowerParBounds=lowerParBoundsPop[[iPop]]
			))
			argsFUpdate <- c( list(Xc, argsFUpdateBlock=argsFUpdateBlock) )
			Xc <- do.call( block$fUpdateBlock, argsFUpdate )
		} #iBlock
		
	}# for iThin0
	
	
	
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
	
	.nPop=2
	Zinit <- ZinitPops <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=4*.nPop, nPop=.nPop)
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
		,argsFLogDen = list(
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
	#mtrace(.updateBlockTwDEMC)
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
	if( 0 == length(block$fLogDen) || !is.function(block$fLogDen)) stop(".checkBlock: entry fLogDen (function of log of unnomalized density) must be specified")
	cpOrig <- block$compPos
	if( is.character(cpOrig) ) block$compPos <- match(cpOrig,parNames)
	if( any(is.na(block$compPos)) ){
		stop(".checkBlock: some of compPos specified in block are not in varNames of initial values.")
	}
	if( 0==length(block$fUpdateBlock) ) block$fUpdateBlock <- .updateBlockTwDEMC
	if( 0==length(block$useMetropolis) ) block$useMetropolis <- TRUE
	if( 0==length(block$intResCompNames) ) block$intResCompNames <- vector("character",0)
	if( 0==length(block$fLogDenScale) ) block$fLogDenScale <- 1
	if( 0==length(block$TFix) ) block$TFix <- vector("numeric",0)
	block$posTFix <- match( names(block$TFix), block$resCompNames )
	
	### block with entries useMetropolis, intResCompNames, fLogDenScale, TFix, posTFix defined
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

.sampleStates <- function(
	### Sample records and initial state from past records.
	Z					##<< population (steps x parms x chains)
	,mZ					##<< number of steps in Z that have been sampled
	,nGenBack			##<< number of samples in recent past, used to generate the proposals
	,nSample			##<< number of proposals to generate
	#,rLogDen
){
	##seealso<< 
	## \code{\link{twDEMCInt}}
	## \code{\link{.generateXPropThin}}
	nChainsPop = dim(Z)[3] 
	nParm = dim(Z)[2]
	##details<<  
	## Random states for chains for difference vector generation are within subsets of populations.
	## This allows simulating independent population of chains.
	## The acceptance rate may differ amonng populations. Hence, the set of previous generations to 
	## randomly select from also differs between poplations.
	# integer array (thinSteps*nChain*4) sample of chains, within populations
	X <- adrop(Z[mZ,,,drop=FALSE],1)		# assume current state X as beginning of the interval
	nStates <- nSample*nChainsPop*3			# need to sample three states for snooker update
		sChains <- 1:nChainsPop 
		sGens <- (mZ-nGenBack+1):mZ 
		rrGenPop <-  sample(sGens, nStates, replace = TRUE)
		rrChainsPop <-  sample(sChains, nStates, replace = TRUE)
		# in order to constrain two dimensions at the same time use the [] subset with an array see ?"["
		rr <- cbind(rrGenPop,rrChainsPop)
		#zLogDen <- array(rLogDen[rr], dim=c(1,nSample*nChainsPop,3), dimnames=list(parms="logDen", steps=NULL, zi=NULL) )
		rrParms <- cbind( rep(rrGenPop,each=nParm), rep(1:nParm, nStates), rep(rrChainsPop,each=nParm) )
		zParms <- array(Z[rrParms], dim=c(nParm,nSample*nChainsPop,3), dimnames=list( parms=rownames(Z), steps=NULL, zi=NULL) )
		# note that parameters are the first dimension		
		#rrParms <- cbind( rep(rrGenPop,nParm), rep(1:nParm, each=nStates), rep(rrChainsPop,nParm) )
		#zParms <- matrix( Z[rrParms], ncol=nParm)
		#z0 <- abind( zParms, zLogDen, along=1)
		
		#append as forth initial state x vector to z
		chainsZ <- rep(sChains,each=nSample)
		Xs <- array(X[,chainsZ], dim=c(nParm,nSample*nChainsPop,1), dimnames=list(parms=rownames(Z), steps=NULL, zi=NULL) )
		#XsLogDen <- array( rLogDen[mZ,chainsZ],dim=c(1,nSample*nChainsPop,1), dimnames=list(parms="logDen", steps=NULL, zi=NULL)  ) 
		#Xs <- abind( Xs, XsLogDen, along=1)	#logDen row never used for X
		zx <- abind( zParms, Xs, along=3)
		
		##details<< 
		## States z2 is  attempted to be distinct.
		## And z3 is attempted to be distinct from x (assuming steps in z are stacked chains collectives of X)  
		iSame <- twWhichColsEqual( zx[,,2], zx[,,1] )
		i <- 0
		while( 0 < length(iSame) && i<5){
			zx[,iSame,2] <- zx[,sample.int(dim(zx)[2],length(iSame)),2]
			iSame <- twWhichColsEqual( zx[,,2], zx[,,1] )
			i<-i+1
		}
		#(tmp <- which( zx[,,2] == zx[,,1], arr.ind=TRUE))
		#when z3 is the same as x, i.e. zx[,,4]
		iSame <- twWhichColsEqual( zx[,,3], zx[,,4] )
		#iSame <- unique(which( (zx[-nrow(zx),,3]==Xs), arr.ind=TRUE )[,2])
		i <- 0
		while( 0 < length(iSame) && i<5){
			zx[,iSame,3] <- zx[,sample.int(dim(zx)[2],length(iSame)),3]
			iSame <- twWhichColsEqual( zx[,,3], zx[,,4] )
			i<-i+1
		}
		#(tmp <- which( zx[,,3] == zx[,,4], arr.ind=TRUE))
		### random states (Nparms+1,(steps*nChainsPop), 4)
		### first dimension is the state vector
		### random states for each step and chains (stacked to be vectorized)
		### chain is last dimensionion in stack (consequtive steps for one chain) in  order to abind across populations
		### zx dim: three random vectors, forth dimension is the initial state x of the chain
		rownames(zx) <- colnames(Z)	# preserve parameter names
		zx
}	

.generateXPropThin <- function(
	### Generate Proposals from given random states.
	zx		##<< output of \code{\link{.sampleStates}} (maybe stackes states)
	,ctrl	##<< list with tuning parameters, see \code{twDEMCBlockInt}
	,nChain	##<< number of chains
){
	##seealso<< 
	## \code{\link{twDEMCInt}}
	## \code{\link{.sampleStates}}
	## \code{\link{.xStepSnooker}}
	## \code{\link{.xStepParallel}}
	
	nParm <- nrow(zx)
	nState <- ncol(zx)
	nSample <- nState/nChain
	#zxl <- lapply( 1:nPop, fzPop ) 
	#zxAndLogDen <- structure( abind(zxl, along=2) , dimnames=list(parms=c(rownames(Z),"logDen"),steps=NULL,zi=NULL))
	#zx <- zxAndLogDen[1:nParm,,,drop=FALSE]
	z <- zx[,,1:3,drop=FALSE]	#three random state vectors per step
	X <- adrop(zx[,,4,drop=FALSE],3)	#initial state vector for step
	#zLogDen <- adrop(zxAndLogDen[,,1:3,drop=FALSE][nParm+1,,,drop=FALSE],1)
	dz <- as.list( structure(dim(z), names=c("parms","steps","zi")) )
	
	res <- matrix( NA_real_, nrow=dz$parms+1, ncol=dz$steps, dimnames=c(list(c(rownames(z),"rExtra")),list(NULL)) )
	boSnooker <- runif(dz$steps)< ctrl$pSnooker
	if( 0 < sum(boSnooker) ){
		res[,boSnooker] <- .xStepSnooker(z[,boSnooker,,drop=FALSE],X[,boSnooker,drop=FALSE],ctrl=ctrl)
	}
	if( 0 < sum(!boSnooker) )
		#res[,!boSnooker] <- .xStepParallel(z[,!boSnooker,,drop=FALSE],zLogDen=zLogDen[!boSnooker,,drop=FALSE],ctrl=ctrl)	
		res[,!boSnooker] <- .xStepParallel(z[,!boSnooker,,drop=FALSE],ctrl=ctrl)	
	
	#second dimension is Nsteps*nChains (we can set nChains as last dimension 
	#array(chainsZ,dim=c(nSample,d$chains))
	resArraySteps <- array(res, dim=c(nParm+1,nSample,nChain), dimnames=list(parms=rownames(res),steps=NULL,chains=NULL) )
	#expected steps as last dimension
	xStepAndExtra <- aperm(resArraySteps,c(1,3,2))
	# numeric array (Npar+1,Nchains,Nsteps): difference vectors in parameter space for steps and chains
	# last row is the extra LogDen associated with snooker update
	
	xStep <- xStepAndExtra[-nrow(xStepAndExtra),,,drop=FALSE]			#dito
	rExtra <- adrop(xStepAndExtra[nrow(xStepAndExtra),,,drop=FALSE],1)			#second dim (columns step within Thinning interval)
	list(xStep=xStep, rExtra=rExtra)
	### List with components \describe{
	### \item{xStep}{numeric array (Npar,Nchain,Nsteps): difference vectors in parameter space}
	### \item{rExtra}{numeric matrix (Npar,Nsteps): some extra LogDensity from snooker update}}
	### Nsteps=ctrl$thin
}

.xStepSnooker <- function(
	### Generates Snooker updates based on given random numbers.
	z, ##<< numeric array (Nparms,(nsteps), 3) of random states, dimnames parms,steps,zi
	X, ##<< current state (Nparms,(nsteps)) corresponding to chain of second dimension in z 
	ctrl
){
	# DE-Snooker update
	##seealso<< 
	## \code{\link{.generateXPropThin}}
	
	nParm <- nrow(z)
	nState <- ncol(z)
	gamma_snooker = runif(nState, min=1.2,max=2.2)
	res <- matrix( as.numeric(NA), nrow=nParm+1, ncol=nState, dimnames=c(list(c(rownames(z),"rExtra")),list(NULL)) )
	#need loop because of inner products
	for( i in 1:nState){
		x_z = X[,i] - z[,i,3]
		D2 = max(1.0e-300, x_z %*% x_z)
		#gamma_snooker =1.7
		projdiff = ((z[,i,1] -z[,i,2]) %*% x_z)/D2  # inner_product of difference with x_z / squared norm x_z
		res[-(nParm+1),i] <- xStepChain <-  (gamma_snooker[i] * projdiff) * x_z
		xPropChain = X[,i] + xStepChain
		x_z = xPropChain - z[,i,3]
		D2prop = max((x_z%*%x_z), 1.0e-30)
		res[(nParm+1),i] <- rExtra <- ctrl$Npar12 * (log(D2prop) - log(D2))   # Npar12  =(nParm - 1)/2  # extra term in logr for accept - reject ratio
	}
	res
} # DE-Snooker update

.xStepParallel <- function(
	### DE-parallel direction update based on given random numbers.
	z, ##<< numeric array (Nparms,(nsteps), 3) of random states, dimnames parms,steps,zi
	#zLogDen,	##<< numeric matrix (nsteps,3): logDen corresponding to the random states z   
	ctrl
){
	##seealso<< 
	## \code{\link{.generateXPropThin}}
	
	nParm <- nrow(z)
	nState <- ncol(z)
	dz <- adrop((z[,,1,drop=FALSE]-z[,,2,drop=FALSE]),3)	#jump vector as the difference between two random states (from z2 towards z1)
	.tmp.f <- function(){
		if( !is.null(ctrl$probUpDir) && (ctrl$probUpDir != 1/2) ){
			iFinites <- which(is.finite(zLogDen[,1]) & is.finite(zLogDen[,2]))
		}
	}
	boGammasF1 <- (runif(nState) < ctrl$pGamma1)
	gamma_par <- matrix( ctrl$F1, nrow=nParm, ncol=nState) # to be able to jump between modes
	gamma_par[,!boGammasF1] <- ctrl$F2 * runif( nParm*sum(!boGammasF1), min=1-ctrl$epsMult, max=1+ctrl$epsMult)    # multiplicative error to be applied to the difference 	
	xStepChain <- if (ctrl$epsAdd ==0) {  # avoid generating normal random variates if possible
			gamma_par * dz 
		} else {
			gamma_par * dz  + rnorm(length(dz),0,ctrl$epsAdd)
		}
	xStepChainAndRExtra <- rbind(xStepChain,rExtra=0)
	### Numeric matrix (Nparms, nsteps): difference vectors in parameter space.
}


.doDEMCSteps <- function(
	### Perform Metropolis steps within next thinning interval.
	X,				##<< numeric matrix: current location of chains rows: parameters, columns: chains
	logDenCompX,	##<< numeric array: result of fLogDen for current location, rows: result components, columns chains  
	#logDenX, 		##<< numeric vector of LogDensity of current position of chains
	xStep, 			##<< array rows: difference vectors in parameter space, cols: chains, zdim: steps
	rExtra,			##<< numeric matrix: row: chain, col: step within thinning interval
	temp,			##<< numeric matrix: row: chain, col: step within thinning interval
	nPop,			##<< the number of populations
	pAccept,		##<< numeric vector: current acceptance rate for each population 
	#argsDEMCStep,	##<< see \code{\link{.doDEMCStep}}
	remoteFunArgs,	
	### see \code{\link[twSnowfall]{sfRemoteWrapper}}
	### Must include entries \itemize{
	### \item{remoteFun=.doDEMCStep}
	### \item{remoteDumpfileBasename} 
	### \item{argsDEMCStep}
	### }
	### Can be a name of a previously exported variable.
	debugSequential=FALSE	##<< see \code{\link[twSnowfall]{sfFArgsApplyDep}}
	,doRecordProposals=FALSE	##<< if TRUE then proposals and results of rLogDen are recorded in result$Y.
){
	# .doDEMCSteps
	##seealso<< 
	## \code{\link{twDEMCInt}}
	## \code{\link{.doDEMCStep}}
	
	##details<< 
	## The step must be the last dimension in all arguments in order to make use of dependence step 
	## in load balanced execution.
	
	#if(!is.numeric(X) | !is.numeric(logDenCompX) | !is.numeric(logDenX) )
	#	stop(".doDEMCSteps: first three arguments must be numeric")
	xStepStacked <- do.call( rbind, lapply(1:(dim(xStep)[3]),function(iStep){t(adrop(xStep[,,iStep,drop=FALSE],3))}) )
	#all chains of first step, all chains of second step, ...
	d <- as.list(structure(dim(xStep),names=c("parms","chains","steps"))) 
	iGenT <- (1:d$steps)
	nCases = d$chains * d$steps
	nChainsPop = d$chains / nPop
	if( !(nCases == length(rExtra)) || !(dim(rExtra)==dim(temp)) )
		stop("number of cases in steps must correspond to length of rExtra and length of temp")
	#iPops <- matrix( 1:d$chains, ncol=nPop)	
	#F_APPLY <- .doDEMCStep	#the function executed on the nodes: one metropolis step
	F_APPLY <- sfRemoteWrapper	#the function executed on the nodes: one metropolis step
	#fArgsSame <- list( remoteFunArgs=as.name("argsUpdateBlocksTwDEMC") )	#exported, includes remoteFun=.doDEMCStep, remoteDumpfileBasename and argsDEMCStep
	F_ARGS <- function(i,prevRes){ 
		iChain0<-((i-1) %% d$chains)
		iPop<-(iChain0 %/% nChainsPop)+1
		args<-c(	
			prevRes[c("x", "logDenCompAcc", "logDenAcc")],
			list( step=xStepStacked[i,,drop=TRUE], rExtra=rExtra[i], temp=temp[i], iPop=iPop, pAccept=pAccept),
			#list( argsDEMCStep=argsDEMCStep )
			list( remoteFunArgs=remoteFunArgs )
		)}	
	#.res0 <- lapply(1:nrow(X),function(row){X[row,]})
	res0 <- lapply(1:d$chains,function(iChain){list(
				x=X[,iChain,drop=TRUE]
				,logDenCompAcc=logDenCompX[,iChain,drop=TRUE]
				,logDenAcc=logDenX[iChain]
			)})
	res <- sfFArgsApplyDep( nCases, F_ARGS, F_APPLY, res0, debugSequential=debugSequential)
	#all chains of first step, all chains of second step, ...
	#iChains <- (d$steps-1)*d$chains+(1:d$chains)	#index of the chains of the last step
	#modify in place, so that dimnames etc are preserved
	endChain0 <- (d$steps-1)*d$chains	#index before last step of all chains
	#s0 <- list(X=X, logDenCompX=logDenCompX, logDenX=logDenX)#save former state
	#boResFLogDenX <- { .nc <- ncol(logDenCompX); ( !is.null(.nc) && (.nc > 0) ) }	#if number of columns > 0
	#logDenCompComp
	for( iChain in 1:d$chains ){
		resChain = res[[endChain0+iChain]]
		X[,iChain] <- resChain$x
		logDenCompX[,iChain] <- resChain$logDenCompAcc
		logDenX[iChain] <- resChain$logDenAcc
	}
	Y <- NULL
	if( doRecordProposals){ 
		parNames <- rownames(X)
		resCompNames <- names(res[[1]]$logDenCompProp)
		if( is.null(resCompNames) ) resCompNames <- paste("logDenComp",1:length(res[[1]]$logDenCompProp),sep="")
		tmp <- c("rLogDen",parNames,"accepted",resCompNames) 
		Y <- array( double(length(tmp)*d$steps*d$chains), dim=c(d=length(tmp),nStep=d$steps,nChain=d$chains),	dimnames=list(comp=tmp,steps=NULL,chains=NULL) )			
		for(iStep in 1:d$steps){
			chain0 <- (iStep-1)*d$chains
			for( iChain in 1:d$chains ){
				i <- chain0+iChain
				Y["accepted",iStep,iChain] <- res[[i]]$accepted
				Y[parNames,iStep,iChain] <- res[[i]]$xProp
				Y[resCompNames,iStep,iChain] <- res[[i]]$logDenCompProp
				Y["rLogDen",iStep,iChain] <- res[[i]]$logDenProp
			}
		}
	}
	acceptedM <- matrix( sapply(res,function(resChain){ resChain$accepted}), byrow=TRUE, ncol=d$chains )	#rows: steps, cols: chains
	accepted <- colSums(acceptedM)
	
	resDo <- list(	X=X, logDenCompX=logDenCompX, logDenX=logDenX, accepted=accepted, Y=Y )
	### list with components \describe{
	### \item{X}{matrix current position, column for each chain}
	### \item{logDenCompX}{matrix: result components of fLogDen current position, column for each chain}
	### \item{logDenX}{vector current logDen of chains}
	### \item{accepted}{numerical vector, number of accepted steps for each chain}
	### \item{Y}{numerical matrix (steps x 2+params+result): accepted, rLogDen, parms, and all fLogDen result components for each proposal }
}

.updateBlocksTwDEMC <- function(
	### perform a Gibbs cycle by calling the update functin for all parameter blocks
){
	# under development directly in twDEMCBlocks	
}

.updateBlockTwDEMC <- function( 
	### Perfrom one DEMC step, function to be called in remote process.
	x				##<< numeric vector: current state
	,argsFUpdateBlock	
	### arguments that do not change between steps: list with components \describe{
	### \item{step}{proposed jump}
	### \item{rExtra}{extra portion in Metropolis decision because of selecting the jump}
	### \item{temp}{global temperature}
	### \item{TFix}{named numeric vector (comp): components with fixed Temperature}
	### \item{posTFix}{integer vector (comp): =match(TFix, compNames): positions of TFix within comp provided for performance reasons}
	### \item{ctrl$useMultiT}{boolean wheter to scale temperature for different data streams}
	### \item{Tprop}{numeric matrix (comp x pops): proportions of temperature for different data streams}
	### \item{fDiscrProp,argsFDiscrProp}{function and additional arguments applied to xProp, e.g. to round it to discrete values}
	### \item{argsFLogDen, fLogDenScale}{additional arguments to fLogDen and scalar factor applied to result of fLogDen}
	### \item{posLogDenInt}{the matching positions of intResCompNames within the the results components that are handled internally}
	### \item{ctrl$DRgamma}{ if !0 and >0 delayed Rejection (DR) (Haario06) is applied by jumping only DRgamma distance along the proposal }
	### \item{upperParBounds}{ named numeric vector, see \code{\link{twDEMCInt}}  }
	### \item{lowerParBounds}{ named numeric vector, see \code{\link{twDEMCInt}}  }
	### \item{fCalcComponentTemp}{ functiont to calculate temperature of result components, (way of transporting calcComponentTemp to remote process) }
	### }
){
	#.updateBlockTwDEMC
	##seealso<< 
	## \code{\link{twDEMCInt}}
	## \code{\link{.updateBlocksTwDEMC}}
	#attach( argsDEMCStep )
	#stop(".doDEMCStep: stop to trace error in remote R invocation.")
	with( argsFUpdateBlock, {
			boResFLogDenX <- (length(posLogDenInt) > 0)
			# LogDensity of accepted state
			#La <- logDenCompAcc	#logDensity	components of accepted state
			La <- fLogDenScale * if(boResFLogDenX){
					logDenCompInt <- logDenCompAcc[posLogDenInt]
					TiInt <- Ti[posLogDenInt]
					do.call( fLogDen, c(list(x, logDenCompInt, TiInt), argsFLogDen) )	# evaluate logDen
				}else
					do.call( fLogDen, c(list(x), argsFLogDen) )	# evaluate logDen
			#assume that all is.finite(logDenCompAcc), make sure in twDEMCInt
			LaExt <- La
			logDenAcc <- sum(La)
			##details<< \describe{\item{Temperature proportions}{
			## If ctrl$useMultiT=TRUE then at high Temperatures, all datastreams are weighted so that 
			## each one has the same influcence.
			## The Temperature of the components with higher LogDensity (less negative)
			## }}
			
			##details<< \describe{\item{Components with fixed temperature}{
			## If ctrl$TFix is given then the given componenents are assigened the given
			## temperature independent of decreasing global temperature.
			## Useful for priors: \code{TFix = c(parms=1)}
			## }}
			
			Ti <- fCalcComponentTemp( temp=temp, TFix=TFix, Tprop=Tprop, useMultiT=ctrl$useMultiT,posTFix=posTFix)
			
			accepted<-FALSE
			xProp = x + step
			
			boOutside <- 
				any( sapply( names(upperParBounds), function(pname){ xProp[pname] > upperParBounds[pname] })) ||
				any( sapply( names(lowerParBounds), function(pname){ xProp[pname] < lowerParBounds[pname] }))
			#if(xProp[1] < 10.8)	recover()
			if( boOutside ){
				# if it is still outside (maybe opposite border) reject step and give -Inf as logDenResult
				logDenProp=logAlpha10=-Inf		#logAlpha10 is log of the initial acceptance ratio for DR step (-Inf no chance of acceptance)
				Lp <- logDenCompAcc	#results for the proposal
				Lp[] <- -Inf
			}else{
				# discrtize proposal
				if( is.function(fDiscrProp)) xProp = do.call(fDiscrProp,xProp,argsFDiscrProp, quote=TRUE)
				Lp <- fLogDenScale * if(boResFLogDenX){
						logDenCompInt <- logDenCompAcc[posLogDenInt]
						TiInt <- Ti[posLogDenInt]
						do.call( fLogDen, c(list(xProp, logDenCompInt, TiInt), argsFLogDen) )	# evaluate logDen
					}else
						do.call( fLogDen, c(list(xProp), argsFLogDen) )	# evaluate logDen
				#take care that the result has always the same sames, even when if fails
				#if( 0==length(names(res)))
				#	stop("encountered result of fLogDen without names")
				#if( !identical(names(logDenCompAcc),names(res)))
				#	stop("encountered result with different names")
				#strip attributes other than names, else twDynamicClusterApplyDep fails with big data chunks
				attributes(Lp) <- list(names=names(Lp))
				logDenProp=-Inf
				#make sure Lp, La have the same order and legnth
				#if( !identical( names(Lp), names(La)) ) stop(".doDEMCStep: logDenCompAcc must contain the same components and the order of result of fLogDen." )
				if( all(is.finite(Lp))){
					logDenProp <- sum(Lp)
					##details<< \describe{\item{internal Metropolis step}{
					## if posLogDenInt is given, then these components of result of fLogDen are handled
					## internally. Hence, for Metropolis step here operates only on other remaining components.
					## }}
					posTExt <- setdiff( seq_along(Ti), posLogDenInt )		#externally handled components
					nExt <- length(posTExt)
					#posTFixExt <- setdiff(posTFix,posLogDenInt)		#externally handled components with fixed temperature
					#posTVarExt <- setdiff(seq_along(Lp), c(posTFix,posLogDenInt))	#externally handled componetns with variable temperature
					#nFixExt <- length(posTFixExt)
					#nVarExt <- length(posTVarExt)
					#nExt <- nFixExt + nVarExt
					
					#Metropolis step
					#logr = (logDenPropExt+rExtra - logDenXExt) / temp
					logrDS10 <- (Lp[posTExt]-La[posTExt])/Ti[posTExt]
					logAlpha10 <- rExtra + sum(logrDS10) 
					accepted <-  is.finite(logAlpha10) && (logAlpha10  > log(runif(1)) )
					if(accepted){
						x <- xProp
						logDenAcc <- logDenProp
						logDenCompAcc <- Lp
					}				
				}else logAlpha10 <- -Inf
			} # end check outside parBounds
			
			if(!accepted && !is.null(ctrl$DRgamma) && (ctrl$DRgamma > 0) && 
				( boOutside ||	(!is.null(argsDEMCStep$minPCompAcceptTempDecr) && (pAccept < 1.2*argsDEMCStep$minPCompAcceptTempDecr)))
				) {
				#----- delayed rejection (DR) step
				# only if across parBoundEdge or acceptance rate drops below 1.2*minAcceptrate
				# repeat all above with delayed rejection (DR) step, only adjust DRfac after calculating Lp
				Lp1 <- Lp
				xProp1 <- xProp
				xProp <- x + ctrl$DRgamma*step
				
				boOutside <- 
					any( sapply( names(upperParBounds), function(pname){ xProp[pname] > upperParBounds[pname] })) ||
					any( sapply( names(lowerParBounds), function(pname){ xProp[pname] < lowerParBounds[pname] }))
				if( !boOutside ){
					if( is.function(fDiscrProp)) xProp = do.call(fDiscrProp,xProp,argsFDiscrProp, quote=TRUE)
					Lp <- fLogDenScale * if(boResFLogDenX){
							logDenCompInt <- logDenCompAcc[posLogDenInt]
							TiInt <- Ti[posLogDenInt]
							do.call( fLogDen, c(list(xProp, logDenCompInt, TiInt), argsFLogDen) )	# evaluate logDen
						}else
							do.call( fLogDen, c(list(xProp), argsFLogDen) )	# evaluate logDen
					#take care that the result has always the same sames, even when if fails
					#if( 0==length(names(res)))
					#	stop("encountered result of fLogDen without names")
					#if( !identical(names(logDenCompAcc),names(res)))
					#	stop("encountered result with different names")
					#strip attributes other than names, else twDynamicClusterApplyDep fails with big data chunks
					attributes(Lp) <- list(names=names(Lp))
					logDenProp=-Inf
					#make sure Lp, La have the same order and legnth
					#if( !identical( names(Lp), names(La)) ) stop(".doDEMCStep: logDenCompAcc must contain the same components and the order of result of fLogDen." )
					if( all(is.finite(Lp))){
						logDenProp <- sum(Lp)
						##details<< \describe{\item{internal Metropolis step}{
						## if posLogDenInt is given, then these components of result of fLogDen are handled
						## internally. Hence, for Metropolis step here operates only on other remaining components.
						## }}
						posTExt <- setdiff( seq_along(Ti), posLogDenInt )		#externally handled components
						nExt <- length(posTExt)
						#posTFixExt <- setdiff(posTFix,posLogDenInt)		#externally handled components with fixed temperature
						#posTVarExt <- setdiff(seq_along(Lp), c(posTFix,posLogDenInt))	#externally handled componetns with variable temperature
						#nFixExt <- length(posTFixExt)
						#nVarExt <- length(posTVarExt)
						#nExt <- nFixExt + nVarExt
						#Metropolis step 
						#logr = (logDenPropExt+rExtra - logDenXExt) / temp
						logrDS20 <- (Lp[posTExt]-La[posTExt])/Ti[posTExt]
						logAlpha20 <- rExtra + sum(logrDS20)
						#---  here correct with first stage DR factor (1-alpha21)/(1-alpha10) with meaning 0:accepted 1:first proposal 2:second proposal
						logrDS21 <- (Lp1[posTExt]-Lp[posTExt])/Ti[posTExt]
						logAlpha21 <- sum(logrDS21)
						logAlpha2 <- suppressWarnings( logAlpha20  +log(1-exp(logAlpha21)) -log(1-exp(logAlpha10)) )	# log and exp may produce NaNs 
						accepted <-  is.finite(logAlpha2) && ( logAlpha2 > log(runif(1)) ) 
						if(accepted){
							x <- xProp
							logDenAcc <- logDenProp
							logDenCompAcc <- Lp
						}				
					}
				} # end !boOutside in DR step 
			}	# end DR step
			
			#will invoke prevRes[c("x", "logDenCompAcc", "logDenAcc")]
			if(!is.numeric(x) | !is.numeric(logDenCompAcc) | !is.numeric(logDenAcc) )
				stop(".doDEMCStep: x, logDenAcc and logDenCompAcc must be numeric")
			
			list(accepted=accepted
				,x=x,logDenCompAcc=La, logDenAcc =logDenAcc		# input to repeated call
				,xProp=xProp,logDenCompProp=Lp, logDenProp=logDenProp
			)
		})
	#detach( argsDEMCStep ); list(accepted=accepted,x=x,logDenCompAcc=logDenCompAcc, logDenAcc =logDenAcc,xProp=xProp,logDenCompProp=res, logDenProp=logDenProp	) 
	### list with components \describe{
	### \item{accepted}{boolean scalar: if step was accepted}
	### \item{x}{numeric vector: current position in parameter space}
	### \item{logDenCompAcc}{numeric vector: result components of fLogDen for current position }
	### \item{logDenAcc}{numeric vector: summed fLogDen for current accepted position}
	### \item{xProp}{numeric vector: proposal}
	### \item{logDenCompProp}{numeric vector: result components of fLogDen for proposal }
	### \item{logDenProp}{numeric vector: summed fLogDen for proposal}
	### }
}

	

