twDEMCBlockInt <- function(
	### Differential Evolution Markov Chain with blocked parameter update
	pops =  list( list( ##<< list of population infos for each population, each info is a list with components
			##describe<< 
			parms		##<< list of matrices (nState x nParm  x nChain) initial states for each population see details and \code{\link{initZtwDEMCNormal}}.
			, nGen = 12	##<< number of generations, i.e. steps to proceed
			, TProp=NULL	##<< numeric vector (nResComp) temperature proportions of result components.
			    ## It can also be given as character vector with names of result components, however, be aware that htis fails if several logDen may return components of the same name
			, X=NULL		##<< numeric matrix (nParm x nChainPop) initial state
			, logDenCompX=NULL 		##<< numeric matrix (nComp x nChain): logDen components of initial state X, see details
			, spaceInd = 1	##<< the space replicate that this population belongs to
			, upperParBounds = numeric(0)   ##<< named numeric vectors: giving upper parameter bounds: lowerBound < par <= upperBound.
			    ## For exploring subspaces of the limiting distribution, see details
			, lowerParBounds = numeric(0)  ##<< similar to upperParBounds: sample > bound
			, splits=numeric(0)		##<< named numeric vector of splitting points, used to remerge divided subspaces 
		)) 
		##end<<  
	, dInfos = list( list(  ##<< named list of used density functions. Each entry is a list with components
			##describe<< 
			fLogDen=NULL                    ##<< \code{function(theta, ...)} calculates a vector of logDensities 
			    ## corresponding to different data streams of parameter vector theta 
			    ## \cr or \code{function(theta, logDenCompX, metropolisStepTemp, ...)}
			,xC=xC						    ##<< numeric vector: current accepted state
			,logDenCompC=logDenCompC	    ##<< numeric vector: result components of fLogDen for current position
			,parUpdateDenC=parUpdateDenC	##<< integer vector: logDensity that recently updated parameter at given index			
			    ## to handle first steps in Multi-step Metropolis decision internally. 
			    ## See details.
			, compPosDen=1:nrow(pops[[1]]$parms)	##<< index or names of the parameter components that are used in this density function 
			, argsFLogDen=list()	                ##<< further arguments passed to fLogDen
			, intResComp=vector("integer",0)        ##<< integer or character vector: indices or names of results components of fLogDen 
			    ## that are used for internal Metropolis decisions
			, fLogDenScale=1                        ##<< scalar multiplied to the result of fLogDen 
			    ## allows using functions of negative LogDensity (-1) or Gaussian misfit function (-1/2) instead of logDensity
			#, TFix = vector("numeric",0) ##<< named numeric vector with Temperature for result components that have fixed temperature
		)) 
		##end<< 
	, blocks = list( list( ##<< list of parameter blocks, each block is a list with entries
			##describe<< 
			dInfoPos=1		##<< name or position to \code{fLogDenInfo}. Several blocks may share the same density but update different parameters
			, compPos=dInfos[[1]]$compPosDen	##<< names or index of the parameter components to be updated
			, fUpdateBlock=updateBlockTwDEMC	##<< function to update the parameters.
			## \cr It must return a list with first three components xC, logDenCompC, and parUpdateDenC 
			## as described in \code{\link{updateBlockTwDEMC}}
			, argsFUpdate=list()	    ##<< further arguments passed to fUpdate
			, requiresUpdatedDen=TRUE	##<< if fUpdateBlock does not depend on current density result, then set this to FALSE and save some calculation time 
			#, useMetropolis=TRUE	##<< TRUE, if jumps for Metropolis proposals should be generated for this block
		))
		##end<<
	, TSpec = matrix(1, ncol=2, nrow=1, dimnames=list(NULL, c("T0","TEnd")) )
		### numeric matrix (nResComp x 2): specifying Initial and End Temperature of each fLogDen result component.
		### If only one row is specified, the Temperature is taken for all result components
	, m0 = numeric(0)			##<< scalar integer: number of samples in initial population, if length==0 will be calculated from number of chains and number of parameters
	, nGen = integer(0)			##<< scalar integer: number of generations, if given overwrites \code{pops[[i]]$nGen}
	, spacePop = 1:nPop			##<< the space replicate that each population belongs to. By default assume only one population per space, overridden by entry in pops 
	, controlTwDEMC = list()	##<< control parameters influencing the update and the convergens of the chains (see details)	
	, debugSequential=FALSE 		##<< if TRUE apply is used instead of sfApply, for easier debugging
	, remoteDumpfileBasename=NULL	##<< the basename of a dumpfile that is created on error on remote process (see example section) 
	, fDiscrProp=NULL				##<< function applied to proposal, e.g. to round proposals to to discrete possible values function(theta,...)
	, argsFDiscrProp=list()			##<< further arguments to fDiscrProp
	, doRecordProposals=FALSE		##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y
	, progressOutput="."			##<< output after each thinning interval
){
	##details<<  
	## This is the central method for applying a Differential Evolution Markov Chain given a function of 
	## an unnormalized probability density.
	## It is invoked usually by (\code{\link{twDEMCBlock.array}} or \code{\link{twDEMCBlock.twDEMCPops}})
	
	##seealso<< 
	## \code{\link{calcDEMCTemp}}
	## \code{\link{logDenGaussian}}
	
	# # \code{\link{.generateXPropThin}}
	# # \code{\link{.doDEMCSteps}}
	# # \code{\link{.doDEMCStep}}
	
	## #details<<
	## # further elements in pop that are calculated are \itemize{
	## #  \item resCompPos: position of the result components in the concatenation across all density functions  
    ## # }
	
	
	##details<< \describe{ \item{recognized control parameters in \code{controlTwDEMC}: }{
	##describe<<
	ctrl = list( 
		thin = 4 	   	##<< thinning interval	 
		,F = 2.38 		##<< related to multiplicative error (F2=F/sqrt(2*Npar), see eps.mult 
		,pSnooker= 0.1	##<< probability of a snooker update (others parallel updates)
		,pGamma1 = 0.1	##<< probability of jumping to state of another chain (different modes)
		,epsMult =0.2	##<< >0 gives d-dimensional variation around gamma. It adds scaled uncorrelated noise to the proposal. Its advantage over eps.add is that its effect scales with the differences of vectors in the population whereas eps.add does not. if the variance of a dimensions is close to 0, eps.mult gives smaller changes. \cr A uniformly distributed error, i.e. F2*runif(1+-epsMult*prop) multiplied to difference vector from parallel update 
		,epsAdd = 0	   	##<< >0 is needed to ensure that all positions in the space can be reached. For targets without gaps, it can set small or even to 0. \cr sd of normally distributed error added to proposal by parallel or snooker update. 
		,pAcceptWindowWidth = 32 ##<< number of generations back over which the acceptance rate is calculated
		,probUpDir=0.5 	##<< probability of direction between two states of increasing Density (increase during burin may accelerate convergence)
		,initialAcceptanceRate=0.25	##<< numeric matrix (nBlock x nChains) initially assumed acceptance rate. Used to calculate the number of generations backwards to sample from
		,DRgamma=0		##<< factor for reducing step length [0..1) in delayed rejection step, 0 means no DR step
		,minPCompAcceptTempDecr=0.15  ##<< if acceptance rate drops below minPCompAcceptTempDecr+0.02 times this level, employ delayed rejection (DR)
		,pIndStep = 1.5 ##<< independent state about after on average about those number of 1.5 accepted steps
		,nPastGen = 10  ##<< factor for determining the number of recent past states to sample during burnin. It is multiplied by the number of parameters. Past generations are calculated by deviding by the number of chains per population
        ,useLoadBalancing=FALSE	##<< if set to TRUE, submits each step separaetely to a free node, else each node gets an entire chain process 
        ,freeMasterNode=FALSE	##<< if set to TRUE, no job is submitted to first node, so that this node can dispatch jobs without waiting,, see \code{\link[twSnowfall]{sfFArgsApplyDep}} 
        #,useMultiT = FALSE	##<< whether to downscale Temperature of result components during burnin
		#moved to block,TFix = vector("numeric",0)		##<< named numeric vector: result components for which temperature shoudl not change		
	)  
	##details<< }}
	ctrl[names(controlTwDEMC)] <- controlTwDEMC
	
	#--  dimensions of pops
	nPop <- length(pops)
	iPops <- 1:nPop
	pops <- if( 0 == length(nGen) ){
			lapply( iPops, function(iPop){ .checkPop( pops[[iPop]], spaceInd=spacePop[iPop] )}) # fill in default values for missing entries 
		}else{
			if( length(nGen)==1) nGen=rep(nGen[1],nPop)
			if( length(nGen) != nPop ) stop("twDEMCBlockInt: lenght of nGen must equal nPop.")
			lapply( iPops, function(iPop){ .checkPop( pops[[iPop]], nGen=nGen[iPop], spaceInd=spacePop[iPop])}) # fill in default values for missing entries 
		}
	ZinitPops <- lapply(pops,"[[","parms")
	parNames <- colnames(ZinitPops[[1]])
	nParm <- ncol(ZinitPops[[1]])
	nChainPop <- dim(ZinitPops[[1]])[3]
	nChain <- nChainPop*nPop 
	popChain <- rep(1:nPop,each=nChainPop)	# population for chain at given index
	chainsPop <- lapply( iPops, function(iPop){ (iPop-1)*nChainPop+(1:nChainPop)}) # chains for given population
	if( 0 == length(m0) ) m0 = calcM0twDEMC(nParm,nChainPop)
	M0Pops <- sapply( ZinitPops, nrow )
	if( any(M0Pops < 0.8*m0) ) warning(paste("twDEMCBlockInt: too few initial cases for populations",paste(which(M0Pops < m0),collapse=",")) )
	nGenPops <- sapply( pops, "[[", "nGen")
	nThinnedGenPops = sapply( nGenPops %/% ctrl$thin, max, 0 )	#number of thinning intervals 
	nGenPops = nThinnedGenPops * ctrl$thin  #number of generations in this batch (do not need to move behind last thinning interval)
	Nz <- M0Pops +(nThinnedGenPops)			   #number of rows in Z after batch
	XPops <- lapply( pops, "[[", "X" )
	XChains <- do.call(cbind,XPops)
	logDenCompXPops <- lapply( pops, "[[", "logDenCompX" ) # see initial states for initializing missing entries
	upperParBoundsPop = lapply( pops, "[[", "upperParBounds" )
	lowerParBoundsPop = lapply( pops, "[[", "lowerParBounds" )
	
	#-- dimensions of fLogDenInfo
	dInfos <- lapply(dInfos, .checkDInfo, parNames=parNames)
	nDen <- length(dInfos)
	iDens <- 1:nDen
	if( is.null(names(dInfos))) names(dInfos) <- paste("dInfo",iDens,sep="")
	tmp <- which(names(dInfos) == "")
	names(dInfos)[tmp] <- paste("dInfo",iDens[tmp],sep="")
	compPosDens <- lapply(dInfos, "[[", "compPosDen") 
	
	#-- dimensions of blocks
	#mtrace(.checkBlock)
	blocks <- lapply( blocks, .checkBlock, dInfos=dInfos, parNames=parNames)
	nBlock <- length(blocks)
	iBlocks <- 1:length(blocks)
	#block <- blocks[[1]]
	dInfoPosBlock <- sapply(blocks, "[[", "dInfoPos") 
	# following refers to position in parNames
	compPosBlock <- lapply( blocks, "[[", "compPos" ) # already transformed to positions in .checkBlock
	# following refers to position in dInfo$compPosDen
	compPosInDenBlock <- lapply( blocks, "[[", "compPosInDen" )
	# make sure that all parameters are updated at least once
	if( 0==length(compPosBlock) || (length(sBlock <- sort(unique(do.call(c, compPosBlock)))) != nParm) || !all.equal( sBlock, 1:nParm, check.attributes = FALSE) ) 
		stop(paste("each parameter (columns of parms) must be updated in at least one block. nParm=",nParm,"updatedVars=",paste(sBlock,collapse=",")) )
	
	#-- initialize further parameters to parallel and snooker steps
	# that depend on pops (nParm, nChainPop)
	ctrl$Npar12  =(nParm - 1)/2  # factor for Metropolis ratio DE Snooker update
	ctrl$F2 = ctrl$F/sqrt(2*nParm)
	ctrl$F1 = 1
	ctrl$gInd <- ctrl$nPastGen*nParm/nChainPop	#number of independent generations to sample from, similarly to M0 (usually ctrl$nPastGen=10)
	# terBraak report that may not sample the distribution, if not using the full past
	# but together with decreasing temperature acceptance rate drops very low
	# hence constrain to the past during burnin
	
	
	#-- evaluate logDen of initial states and get names of the result components (dInfo$resCompNames)
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
					#all( is.finite(logDenCompXPops[[iPop]][,1]))
					#chainIsAllFinite <- apply( logDenCompXPops[[iPop]], 2, all, is.finite )
					chainIsAllFinite <- apply( logDenCompXPops[[iPop]], 2, function(cChain){ all(is.finite(cChain))} )
					iChainInPop=which(!chainIsAllFinite)
					if( 0 != length(iChainInPop) )
						cbind(iPop=iPop,iChainInPop=iChainInPop)
					else
						cbind(iPop=iPop, iChainInPop=NA_integer_ )[FALSE,]
				})))
	if( 0 < nrow(missingPopInfo) ){
		# evaluate logDen for those
		missingPopInfo <- cbind(missingPopInfo, iChain=(missingPopInfo[,"iPop"]-1)*nChainPop + missingPopInfo[,"iChainInPop"])
		XChainsMissing <- XChains[,missingPopInfo[,"iChain"] ,drop=FALSE]	# initial states of the missings
		logDenCompDen <- vector("list", nDen )		
		for( iDen in iDens ){
			dInfo = dInfos[[iDen]]
			#tmp <- dInfo$fLogDen; mtrace(tmp)
			#tmp2 <- do.call( tmp, c(list(XChainsMissing[,1]),dInfo$argsFLogDen))
			#mtrace(twCalcLogDenPar)
			.resLogDenPar <- twCalcLogDenPar(
				fLogDen=dInfo$fLogDen
				, xProp=t(XChainsMissing[,,drop=FALSE])
				, logDenCompX=numeric(0)
				, argsFLogDen=dInfo$argsFLogDen
				, fLogDenScale=dInfo$fLogDenScale	
				#no internal decisions, intResCompNames=dInfo$intResCompNames	
				, debugSequential=debugSequential
				, remoteDumpfileBasename=remoteDumpfileBasename 
			)
			if( any(!is.finite(.resLogDenPar$logDen)) )
				stop(paste("twDEMCBlockInt: non-finite logDensity of starting value ",sep=""))
			dInfos[[ iDen ]]$resCompNames <- rcNames <- .getResFLogDenNames( t(.resLogDenPar$logDenComp) ) 
			#dInfos[[ iDen ]]$resCompNamesUnique <- paste(rcNames,iDen,sep="_")  
			#dInfos[[ iDen ]]$resCompNamesUnique <- rcNames  	# keep original names
			logDenCompDen[[iDen]] <- .resLogDenPar$logDenComp
		} #iDen
		logDenCompXMiss <- do.call( cbind, logDenCompDen )
		for( jPop in seq_along(iMissingPops) ){
			iiPop <- iMissingPops[jPop]
			iInfo <- which(  missingPopInfo[,"iPop"] == iiPop )
			logDenCompXPops[[ iMissingPops[jPop] ]] <- 
				t(logDenCompXMiss[iInfo,])
		}
		for( jPop in seq_along(iNonMissingPops) ){
			iiPop <- iNonMissingPops[jPop]
			iInfo <- which(  missingPopInfo[,"iPop"] == iiPop )
			logDenCompXPops[[iiPop]][,missingPopInfo$iChainInPop[iInfo] ] <-
				t(logDenCompXMiss[iInfo,])
		}
	}else{ # any isMissing
		# no Missings, evaluate fLogDen to obtain result component names
		for( iDen in iDens ){
			dInfo = dInfos[[iDen]]
			.resFLogDen <- do.call( dInfo$fLogDen, c(list(XChains[,1]),dInfo$argsFLogDen))
			dInfos[[ iDen ]]$resCompNames <- rcNames <- .getResFLogDenNames( .resFLogDen ) 
			#dInfos[[ iDen ]]$resCompNamesUnique <- rcNames #paste(rcNames,iDen,sep="_")  
		} #iDen
	}
	logDenCompXChains <- do.call( cbind, logDenCompXPops )
	
	#-- make unique resCompNames and calculate positions in resCompNames of results of different blocks
	resCompNamesOrig <- lapply( dInfos, "[[", "resCompNames" )
	resCompNamesU <- resCompNamesOrig #lapply( dInfos, "[[", "resCompNamesUnique" )
	resCompNamesUFlat <- as.vector(do.call(c, resCompNamesU))	# concatenate names
	nResComp <- length(resCompNamesUFlat)
	resCompPosPops <- { # list: for each dInfo: position of resCompNames in resCompNamesUFlat
		tmp.length <- sapply(resCompNamesU, length)
		tmp.end <- cumsum(tmp.length)
		tmp.i <- vector("list",nDen)
		for( iDen in iDens ){ 
			tmp.i[[iDen]] <- dInfos[[iDen]]$resCompPos <- tmp.end[iDen]+1-(tmp.length[iDen]:1) 
		}
		tmp.i
	}
	# update the names in intial conditons to hold unique names
	rownames(logDenCompXChains) <- resCompNamesUFlat
	
	#-- get indices of internal components and TFix (dInfo$ intResCompPosWithin, intResCompPos, posTFix) 
	# depends on proper resCompNames (inferred from invoking fLogDen with initial states)
	intResCompDen <- vector("list",nDen)
	#.TFixL <- vector("list", nDen )
	#.posTFixL <- vector("list", nDen )
	for( iDen in iDens ){
		dInfo = dInfos[[iDen]]
		#.TFixL[[iDen]] <- dInfo$TFix
		#.posTFixL[[iDen]] <- dInfos[[iDen]]$posTFix <- dInfo$resCompPos[ match( names(dInfo$TFix), dInfo$resCompNames ) ] # positions in concatenated result components
		intResCompI <- dInfo$intResComp
		iPosW <- if( is.character(intResCompI) ){
				iPosW <- match(intResCompI, dInfo$resCompNames )				
				if( any(is.na(iPosW)))
					stop(paste("not all names of intResCompNames (",paste(intResCompI,collapse=","),") in return of fLogDen: ",paste(dInfo$resCompNames,collapse=",")))
				iPosW
			} else intResCompI
		dInfos[[iDen]]$intResCompPosWithin <- iPosW
		intResCompDen[[iDen]] <- dInfos[[iDen]]$intResCompPos <- dInfo$resCompPos[iPosW] # position in concatenation of all result components
	} #iDen
	#TFixAll <- do.call(c, .TFixL )
	#posTFixAll <- do.call(c, .posTFixL )	# position in concatenation of all resComp
	
	#-- Temperature of result components, i.e. data streams
	if( nrow(TSpec) == 1 ){
		TSpec <- do.call( rbind, lapply(1:nResComp,function(i){TSpec[1,]}) )
		rownames(TSpec) <- resCompNamesUFlat
	}
	tempResCompPops <- vector("list", nPop )
	for( iPop in iPops){
		.tempResL <- lapply( 1:nResComp, function(iResComp){
			.tempG <- pmax(1,c( TSpec[iResComp,"T0"], 
				calcDEMCTemp( TSpec[iResComp,"T0"], TSpec[iResComp,"TEnd"], nGenPops[iPop] )))
		})
		tempResCompPops[[iPop]] <- do.call( cbind, .tempResL)
		colnames(tempResCompPops[[iPop]]) <- resCompNamesUFlat
	}
	
	#-- setup acceptance rate recording
	# number of requried independent generations to choose from (10*d independent states: TerBraak06 TerBraak08) devided by number of chains
	if( length(ctrl$initialAcceptanceRate) == 1)
		ctrl$initialAcceptanceRate <- matrix(ctrl$initialAcceptanceRate, nrow=nBlock, ncol=nChain, dimnames=list(blocks=NULL, chains=NULL) )
	if( !is.matrix(ctrl$initialAcceptanceRate) || dim(ctrl$initialAcceptanceRate) != c(nBlock,nChain))
		stop("ctrl$initialAcceptanceRate must be matrix of dim nBlock x nChain")
	pAcceptPops <- popMeansTwDEMC(ctrl$initialAcceptanceRate, nPop) # current acceptance rate of population (block x pop)
	pAcceptPops[] <- pmax(1/ctrl$thin,pAcceptPops ) # lower bound
	#nGenBack <- pmin(mZ,ceiling(ctrl$rIndToAcc * ctrl$gInd * pmax(1,ctrl$pIndStep/(ar * ctrl$thin))))	#number of genrations to select states from for each population, ctrl$gInd is multiplied by number of rows for one step depending on acceptance rate and thinning but at least one
	pAcceptMinPops <- apply(pAcceptPops,2,min)	# minimum of acceptance rate across blocks for each population
	nGenBackPops <- pmin(M0Pops,ceiling(ctrl$gInd * pmax(1,ctrl$pIndStep/(pAcceptMinPops * ctrl$thin))))	#number of genrations to select states from for each population, ctrl$gInd is multiplied by number of rows for one step depending on acceptance rate and thinning but at least one  
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
	# (stepslots x blocks x chains)
	#acceptWindow <- matrix( rep(ctrl$thin*ar, each=2*aThinWinW*nChainPop), nrow=2*aThinWinW, ncol=d$chains )	#record number of accepted steps in thinning interval
	# (width x block x chain )
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
	# several resLogDen per chain (dInfos are concatenated)
	resLogDen <- lapply(iPops,function(iPop){ res <- array( NA_real_
				, dim=c(1+nThinnedGenPops[iPop], nResComp, nChainPop)
				, dimnames=list(steps=NULL, resComp=resCompNamesUFlat, chains=NULL))
			res[1,,] <- logDenCompXPops[[iPop]]
			res
		}) 
	# one rLogDen per chain and per dInfo
	logDen <- lapply(iPops, function(iPop){	res <- array(NA_real_
				, dim=c(1+nThinnedGenPops[iPop], nDen, nChainPop) 
				, dimnames=list(steps=NULL, den=names(dInfos), chains=NULL) )
			for( iChainPop in 1:nChainPop)
				for( iDen in iDens )
					res[ 1,iDen,iChainPop] <- sum(logDenCompXPops[[iPop]][resCompPosPops[[iDen]],iChainPop ])
			res
		})
	# one pAcceptance indication per chain and per block
	pAccept <- lapply(iPops, function(iPop){	res <- array(NA_real_
				, dim=c(1+nThinnedGenPops[iPop], nBlock, nChainPop) 
				, dimnames=list(steps=NULL, block=names(blocks), chains=NULL) )
		})
	for( iPop in iPops ){
		pAccept[[iPop]][1,,] <- ctrl$initialAcceptanceRate[,chainsPop[[iPop]] ]
	}
	# record of proposals and fLogDen results, rows c(accepted, parNames, resCompNames, rLogDen)
	# if doRecordProposals is FALSE record only the thinning intervals for the last 128 generations
	# +1 for the initial state
	nThinLastPops <- if(doRecordProposals) nThinnedGenPops else pmin(nThinnedGenPops, ceiling(128/ctrl$thin)) #number of thinning steps to record
	nThinOmitRecordPops = nThinnedGenPops-nThinLastPops	#the Thinning intervals with no recording of outputs
	nGenOmitRecordPops = nThinOmitRecordPops*ctrl$thin
	nGenYPops <- nThinLastPops*ctrl$thin #+1?
	YPops <- lapply(iPops, function(iPop){ res <- array( NA_real_ 
				, dim=c( nGenYPops[iPop], nParm+nBlock+nResComp, nChainPop)
				, dimnames=list(steps=NULL, vars=c(rownames(XPops[[iPop]]),paste("accepted",iBlocks,sep=""),colnames(resLogDen[[iPop]])), chains=NULL ) )
			# already initialized first row in nGen cycle
			#if( nGenYPops[iPop] == nGenPops[iPop] ) # record first state as accepted
			#	res[1,,] <- rbind( XPops[[iPop]], matrix(TRUE, nrow=nBlock, ncol=nChainPop), logDenCompXPops[[iPop]] )
			res
		})	# records proposals, indication of acceptance per Block for each Chain
	# here accepted is 0: not accepted, 1: full step accepted, between DR step accepted at DRGamma		 
	
	#-- arguments to fUpdate argument argsFUpdate that do not change within thinning interval
	argsFUpdateDefault = list(
		ctrl=ctrl
		#,fCalcComponentTemp=calcComponentTemp 
		,fDiscrProp=fDiscrProp
		,argsFDiscrProp=argsFDiscrProp
	)
	# construct an list for each block that includes compPos, and intResCompPos
	argsFUpdateBlocks <- lapply( blocks, function(block){ 
			dInfo <- dInfos[[block$dInfoPos]]
			c( argsFUpdateDefault, dInfo, block
			#, list(TPropPops=TPropPops[dInfo$compPosDen,]) 
			)
		})
	remoteArgsFUpdateBlocksTwDEMC <- list( 
		remoteFun=.updateBlocksTwDEMC,	# this function is called remotely and invokes the update function for each block 
		argsUpdateBlocksTwDEMC=list(
			argsFUpdateBlocks = argsFUpdateBlocks
			,nBlock=nBlock #, iBlocks=iBlocks
			,upperParBoundsPop=upperParBoundsPop
			,lowerParBoundsPop=lowerParBoundsPop
		#,temp=temp
		)
	)
	remoteArgsFUpdateBlocksTwDEMC$remoteDumpfileBasename<-remoteDumpfileBasename	#if null delete
	if( !debugSequential & sfParallel() )
		sfExport("remoteArgsFUpdateBlocksTwDEMC")	#evaluated in remote process, only passed one time
	#else( assign("remoteArgsFUpdateBlocksTwDEMC", remoteArgsFUpdateBlocksTwDEMC, pos=1))	#export to current 
	tmp.remoteFunArgs <- if( !debugSequential && sfParallel() ) as.name("remoteArgsFUpdateBlocksTwDEMC") else remoteArgsFUpdateBlocksTwDEMC 	# eval.parent fails for remoteArgsFUpdateBlocksTwDEMC within sequential function
	
	# arguments to demc that change between thinning intervals
	# substract logDen components of logDenCompX from logDenXExt
	chainState <- list(
		X = XChains
		,logDenCompX = logDenCompXChains
		,parUpdateDen = array(TRUE, dim=c(nDen,nParm,nChain), dimnames=list(dens=NULL,parms=parNames,chains=NULL))  ##<< for each parameter/density combination: is the density up to date
	#,logDenX = logDenX
	)
	parUpdateDenLast = chainState$parUpdateDen # to store the last parUpdateDen
	
	maxNThinned <- max(nThinnedGenPops)
	isSamePopLength = all( nThinnedGenPops == maxNThinned)
	isPops <- 1:nPop	# populations that take part in this step (some may need less steps)
	isChains <- 1:nChain  # chains that take part in this step
	#iThin0 <- maxNThinned-1
	
	#--- main cycle of updates (each ctrl$thin generations) 
	for( iThin0 in (0:(maxNThinned-1)) ){
		if( length(progressOutput)!=0 && nchar(progressOutput)!=0 ){ cat(progressOutput) }
		iGen = iThin0*ctrl$thin+(1:ctrl$thin)	# the generations within this thinning step
		mZPops <- M0Pops + iThin0
		iPopsOut <- which(nThinnedGenPops == iThin0)	# those pops drop out
		if( 0 != length(iPopsOut) ){
			iChainsOut <- do.call( c, lapply( iPopsOut, function(iPop){ chainsPop[[iPop]] }) )
			iiPopsOut <- which(isPops %in% iPopsOut ) # index in set of current populations
			iiChainsOut <- rep((iiPopsOut-1)*nChainPop, each=nChainPop) + (1:nChainPop)
			isPops <- isPops[-iiPopsOut]
			isChains <- isChains[-iiChainsOut]
			# recored the last parUpdateDen for dropouts
			parUpdateDenLast[,,iChainsOut] <- chainState$parUpdateDen[,,iiChainsOut]
			# adapt chainState that refer to current populations only
			chainState <- within(chainState,{
					X <- X[,-iiChainsOut ,drop=FALSE]
					logDenCompX <- logDenCompX[,-iiChainsOut ,drop=FALSE]
					parUpdateDen <- parUpdateDen[,,-iiChainsOut ,drop=FALSE]
				})
			#acceptance rates and nGenBack refer to all chains
		}	
		
		# sample random vectors (note here parameters in rows) 
		# calculate proposed steps (differences not destinations) within next thinning interval
		zxl <- lapply( iPops[isPops], function(iPop){
				.sampleStates(Z=ZPops[[iPop]],mZ=mZPops[iPop],nGenBack=nGenBackPops[[iPop]], nSample=ctrl$thin )				
			})	
		zx <- abind( zxl, along= 2 )	# combine chains (and steps within each chain) of all populations
		genPropRes <- .generateXPropThin(zx,ctrl=ctrl, nChain=length(isPops)*nChainPop)
		
		# check if proposals should be recorded
		# for the first time
		iPopsRecord <- which( (1:nPop %in% isPops) & (iThin0 == nThinOmitRecordPops))
		for( iiPop in iPopsRecord ){
			# record the initial state of the proposals record
			iisPop <- match(iiPop,isPops)	# index in iPop 
			YPops[[iiPop]][1,seq_along(parNames),] <- chainState$X[ ,chainsPop[[ iisPop ]] ]
			YPops[[iiPop]][1,nParm+iBlocks,] <- 1	#TRUE
			YPops[[iiPop]][1,nParm+nBlock+seq_along(resCompNamesUFlat),] <- chainState$logDenCompX[ ,chainsPop[[ iisPop ]] ]
		}
		# after the proposal has been made
		isRecordProposalsPop = ((iThin0 >= nThinOmitRecordPops))[isPops]	# only record the last proposals after nThinOmitRecord
		
		# numeric matrix (nGenThin x nPop)
		#tempGlobalThinStepsL <- lapply( iPops[isPops], function(iPop){ tempGlobalPops[[iPop]][ 1+iGen ] })
		#tempGlobalThinSteps <- structure( do.call( cbind,tempGlobalThinStepsL ), dimnames=list(steps=NULL,pops=NULL) )	
		# numeric array (nGenThin x nResComp x nPop)		
		tempDenCompThinStepsL <- lapply( iPops[isPops], function(iPop){ tempResCompPops[[iPop]][1+iGen, ,drop=FALSE] })
		tempDenCompThinSteps <- structure( abind( tempDenCompThinStepsL, along=3 ), dimnames=list(steps=NULL,resComp=resCompNamesUFlat, pops=NULL) )	
		
		# here may use code in .tmp.f.testStep
		
		#-- do the steps of next thinning interval in load balanced way
        fUpdateInterval <- if( ctrl$useLoadBalancing ) .updateIntervalTwDEMCPar else  .updateIntervalTwDEMCParChains
        # fUpdateInterval <- .updateIntervalTwDEMCParChains
        resUpdate <- fUpdateInterval( X=chainState$X, logDenCompX=chainState$logDenCompX, parUpdateDen=chainState$parUpdateDen
			,xStep=genPropRes$xStep, rExtra=genPropRes$rExtra
			#, tempGlobalSteps=tempGlobalThinSteps
			, tempDenCompSteps=tempDenCompThinSteps		
			, nsPop=length(isPops), pAcceptPar=pAcceptPops[,isPops ,drop=FALSE]
			,remoteFunArgs=remoteArgsFUpdateBlocksTwDEMC
			,debugSequential=debugSequential
			,isRecordProposalsPop=isRecordProposalsPop 
			,isPops=isPops
			,freeMasterNode=ctrl$freeMasterNode
		)
		chainState <- resUpdate[ names(chainState) ]
		
		#-- calculate accepteance rate
		#row in acceptance Window to record acceptance, if exceeds window, copy second part to first (rewind)
		acceptPos0 <- aThinWinW + (iThin0 %% aThinWinW)
		if( acceptPos0 == aThinWinW ){	#interval exeeded or new, copy second half to first half
			acceptWindow[ 1:aThinWinW,, ] <- acceptWindow[ aThinWinW+(1:aThinWinW),, ]
			acceptWindow[ aThinWinW+(1:aThinWinW),, ] <- NA			
		}
		acceptWindow[ acceptPos0+1,, isChains] <- resUpdate$accepted
		curAcceptRows <- (acceptPos0+1)-(min(aThinWinW-1,max(iThin0,8)):0)	# at the start average over at least 4 slots,  
		#acceptWindow[curAcceptRows,]
		pAcceptChains <- apply(acceptWindow[curAcceptRows,, ,drop=FALSE],c(2,3),sum) / (length(curAcceptRows)*ctrl$thin)
		#if( !all(is.finite(pAcceptChains[,isChains ,drop=FALSE])))
		#	stop("twDEMCBlocks: encountered non-finite pAcceptChains")
		pAcceptPops <- popMeansTwDEMC(pAcceptChains, nPop) # current acceptance rate of population (block x pop)
		pAcceptPops[] <- pmax(1/ctrl$thin,pAcceptPops ) # lower bound
		pAcceptMinPops <- apply(pAcceptPops,2,min)	# minimum of acceptance rate across blocks for each population 
		nGenBackPops <- pmin(M0Pops,ceiling(ctrl$gInd * pmax(1,ctrl$pIndStep/(pAcceptMinPops * ctrl$thin))))	#number of genrations to select states from for each population, ctrl$gInd is multiplied by number of rows for one step depending on acceptance rate and thinning but at least one  
		
		#-- record current thinning step for each population
		for( iiPop in seq_along(isPops) ){
			iPop <- isPops[iiPop]		# translate index in current populations to all populations 
			mZ <- mZPops[ iPop ]
			ipChains <- chainsPop[[iiPop]]	# chains within current chains
			igChains <- chainsPop[[iPop]]	# chains within all chains
			ZPops[[iPop]][mZ+1,,] <- chainState$X[,ipChains]
			resLogDen[[iPop]][iThin0+2,,] <- chainState$logDenCompX[,ipChains] # first state is initial state
			pAccept[[iPop]][iThin0+2,,] <- pAcceptChains[,igChains]
#if( iPop==13 && any(chainState$X["a",ipChains] < -0.6) ) recover()	
		}
		
		#-- record poprosals
		iYSteps0 <- (iThin0 - nThinOmitRecordPops)*ctrl$thin 
		for( iiPop in which(isRecordProposalsPop) ){
			# recordProposals refers to current populations (some may have dropped out already)
			iPop <- isPops[iiPop]
			iiChains = chainsPop[[iiPop]]
			YPops[[iPop]][iYSteps0[iPop]+(1:ctrl$thin),,] <- resUpdate$Y[,,iiChains]
		}
	}# for iThin0
	if( length(progressOutput)!=0 && nchar(progressOutput)!=0 ){ cat("\n") }
	# recored the last parUpdateDen for dropouts
	iChainsOut <- do.call( c, lapply( isPops, function(iPop){ chainsPop[[iPop]] }) )
	parUpdateDenLast[,,iChainsOut] <- chainState$parUpdateDen
	
	#-- calculate the log-Density of last states that are not up to date, because they might be used for reinitialization
	for( iDen in iDens ){
		dInfo <- dInfos[[iDen]]
		cd <- dInfo$compPosDen # components used by the current density function
		udi <- adrop(parUpdateDenLast[iDen,, ,drop=FALSE],1)   # (iDen= x par x chain
		iChains <- which( apply(udi,2,function(ud){ !all(ud) }) )  
		# TODO parallelize those calculations
		for( iChain in iChains ){
			iPop <- popChain[iChain]
			iChainInPop <- ((iChain-1) %% nChainPop) +1
			Xc <- ZPops[[iPop]][(M0Pops+nThinnedGenPops)[iPop],cd,iChainInPop ]
			resCompDenC <- dInfo$fLogDenScale * do.call( dInfo$fLogDen, c(list(Xc), dInfo$argsFLogDen) )	# evaluate logDen  
			resLogDen[[iPop]][ 1+nThinnedGenPops[iPop],dInfo$resCompPos,iChainInPop] <- resCompDenC
		} # end need update
	} # for iDen
	
	#-- sum over resLogDen components to yield logDens for each given density 
	# (might be slightly off because not up to date, but work for inspecting the trend)
	for( iDen in iDens ){
		dInfo <- dInfos[[iDen]]
		rcp <- dInfo$resCompPos # components used by the current density function
		if( length(rcp) == 1){
			for( iPop in iPops ){
				logDen[[iPop]][,iDen,] <- resLogDen[[iPop]][,rcp,]
			}
		}else{
			for( iPop in iPops ){
				resLogDenI <- resLogDen[[iPop]][,rcp, ,drop=FALSE]
				logDen[[iPop]][,iDen,] <- apply(resLogDenI,c(1,3),sum)
			}
		} # length(resComp) == 1 
	} # iDen
	
	#-- return
	spacesPop <- sapply(pops,"[[","spaceInd")
	##value<< a list with entries of populations, each entry is a list
	res <- list(
		thin=ctrl$thin		##<< thinning interval that has been used
		,dInfos=dInfos		##<< list of information on densities (argument \code{dInfos})
		,blocks=blocks		##<< list of information on blocks (argument \code{blocks})
		,YPos = list(		##<< list of column positions in Y, a list with entries \describe{
			##describe<<
			accepted = nParm+iBlocks	##<< integer vector of positions of acceptance indication of block at given index
			,resLogDen0 = nParm+nBlock 	##<< integer scalar: postion before first column of results of fLogDen
			)
			##end<<
		,pops = lapply( iPops, function(iPop){ list(  ##<< info on each population. A list with entries: \describe{   
					##describe<<
					upperParBounds = upperParBoundsPop[[iPop]]	##<< upper parameter bounds for sampling
					,lowerParBounds = lowerParBoundsPop[[iPop]] ##<< lower parameter bounds for sampling
					,splits=pops[[iPop]]$splits		##<< named numeric vector: splitting history
					,spaceInd = spacesPop[iPop]		##<< the space replicate that the population belongs to
					,parms = ZPops[[iPop]][ M0Pops[iPop]:nrow(ZPops[[iPop]]),, ,drop=FALSE]	##<< numeric array (steps x parms x chains): collected states, including the initial states
					,temp = tempResCompPops[[iPop]][seq(1,nGenPops[iPop]+1,by=ctrl$thin), ,drop=FALSE] ##<< numeric array (nSample+1, nResComp): temperature, i.e. cost reduction factor in each step for each datastream
					,pAccept= pAccept[[iPop]]	##<< acceptance rate of chains (nStep x nChainPop)
					,resLogDen = resLogDen[[iPop]]	##<< numeric array (steps x resComps x chains): results components of fLogDen of blocks  
					,logDen = logDen[[iPop]]	##<< numberic array (steps x iDen x chains): results summed over blocks
					,Y = YPops[[iPop]]
				)}) 
				##end<<
	) 
	##end<<
	#names(res$pops) <- names(pops)
	class(res) <- c( class(res), "twDEMCPops" )	#monte carlo proposal list
	res
	
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
		)
	)
	#tmp <- .checkPop(pops[[1]])
	
	# both blocks compare the same model against the same observations
	blockDefault <- list(
		fLogDen=logDenGaussian
		,resCompNames=c("obs","parms")
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
	#mtrace(updateBlockTwDEMC)
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
	
	# tracing error in remote session:
	# Provide argument \code{remoteDumpfileBasename="dumpFile"}
	# Then a dumpfile is created on remote error by \code{\link{sfRemoteWrapper}}.
	# In order to trace the density function, the following can be done
	if( FALSE ){
		.remoteDumpfileBasename="dumpfile"
		.remoteDumpfile <- paste(.remoteDumpfileBasename,".rda",sep="")
		load(.remoteDumpfile)
		debugger(get(.remoteDumpfileBasename))
		# choose last step (18)
		require(debug)
		fDen <- remoteFunArgs$argsUpdateBlocksTwDEMC$argsFUpdateBlocks[[1]]$fLogDen
		mtrace(fDen); remoteFunArgs$argsUpdateBlocksTwDEMC$argsFUpdateBlocks[[1]]$fLogDen <- fDen
		do.call(remoteFun, c(remoteFunArgs, list(...)))	
		# go()
		# qqq()
	}
	
	detach()
	
}

.tmp.f.testStep <- function(){
	# test one step for one chain
	iStep <- 1
	iChain <- 1
	# those are arguments, generated updateBlocksTwDEMCPar
	iPop <- popChain[iChain]
	step <- genPropRes$xStep[,iChain,iStep]
	rExtra <- genPropRes$rExtra[iChain,iStep]
	Xc <- chainState$X[,iChain]
	logDenCompC <- chainState$logDenCompX[,iChain]
	parUpdateDenC <- chainState$parUpdateDen[,,iChain]
	tempC <- tempThinSteps[iChain,iStep]
	.updateBlocksTwDEMC(iPop = iPop,step=step,rExtra=rExtra,Xc=Xc,logDenCompC=logDenCompC,parUpdateDenC=parUpdateDenC,pAccept=pAccept, temp=tempC, argsUpdateBlocksTwDEMC=remoteArgsFUpdateBlocksTwDEMC$argsUpdateBlocksTwDEMC )			
}

.tmpf. <- function(){
	# before (code here) calculated temperature for each chain, but actually need it only by population
	tempThinStepsL <- lapply( iPops, function(iPop){
			matrix( temp[[iPop]][ 1+iGen ], nrow=length(iGen), ncol=nChainPop )
		})
	tempThinSteps <- t(do.call( cbind, tempThinStepsL))	# steps need to be in last dimension (so chains in first)
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
	,nGen=numeric(0)		##<< number of generations to take: overwrites settings given in pop
	,spaceInd=numeric(0)	##<< the space replicate that this population belongs to, overwritten by setting in pop 
){
	##details<<
	## If several entries are null, they are set to default values
	if( 0==length(pop$parms) || !is.numeric(pop$parms) || length(dim(pop$parms)) != 3 )
		stop(".checkPop: entry parms must be a numeric 3D array")
	if( 0 != length(nGen) )	pop$nGen <- nGen 
	if( (1 != length(pop$nGen)) )	
		stop(".checkPop: scalar integer nGen must be given as entry nGen in population or in invocation of twDEMCBlockInt.")
	if( 0 == length(pop$spaceInd) ) pop$spaceInd = spaceInd
	if( (1 != length(pop$spaceInd)) )	
		stop(".checkPop: scalar integer spaceInd must be given as entry nGen in population or in invocation of twDEMCBlockInt.")
	# if( 0==length(pop$T0) ) pop$T0 <- 1
	# if( 0==length(pop$TEnd) ) pop$TEnd <- 1
	if( 0==length(pop$X) ){
		pop$X <- adrop( pop$parms[nrow(pop$parms),, ,drop=FALSE],1)
	}else{
		if( !all.equal( rownames(pop$X), colnames(pop$parms), check.attributes = FALSE) )
			stop(".checkPop: rownames of X and colnames of parms must match")
	}
	if( 0 != length(pop$upperParBounds) && any(is.na(pop$upperParBounds)))
		stop(".checkPop: encountered NA in upperParBounds")
	if( 0 != length(pop$lowerParBounds) && any(is.na(pop$lowerParBounds)))
		stop(".checkPop: encountered NA in lowerParBounds")
	
	
	### pop with entries nGen, and X defined 
	pop	
}

.checkDInfo <- function(
	### filling default values in logDenInfos description
	dInfo	##<< list to be checked for entries
	,parNames	##<< variable names
){
	if( 0 == length(dInfo$fLogDen) || !is.function(dInfo$fLogDen)) stop(".checkLogDenInfo: entry fLogDen (function of log of unnomalized density) must be specified")
	##details<< if compPosDen is not given, then assume all parameter dimensions are required	
	if( 0 == length(dInfo$compPosDen) ) dInfo$compPosDen <- seq_along(parNames)
	cpOrig <- dInfo$compPosDen
	if( is.character(cpOrig) ) dInfo$compPosDen <- match(cpOrig,parNames)
	if( any(is.na(dInfo$compPosDen)) ){
		stop(".checkLogDenInfo: some of compPosDen specified in fLogDenInfos are not in varNames of initial values.")
	}
	if( 0==length(dInfo$intResComp) ) dInfo$intResComp <- vector("character",0)
	if( 0==length(dInfo$fLogDenScale) ) dInfo$fLogDenScale <- 1
	#if( 0==length(dInfo$TFix) ) dInfo$TFix <- vector("numeric",0)
	if( 0!=length(dInfo$TFix) ) warning("Usage of TFix is deprecated. Use TSpec with T0 and TEnd of one for given components.")
	if( length(i <- which(!(names( dInfo) %in% c("fLogDen","xC","logDenCompC","parUpdateDenC","compPosDen","intResComp","argsFLogDen","fLogDenScale","intermediate"
                                                ,"resCompNames","resCompPos","intResCompPos","intResCompPosWithin")))) )
        stop(paste(".checkDInfo: unknown entries in dInfo:",paste(names(dInfo)[i],collapse=",") ))
	### argument \code{fLogDenInfo} with entries fLogDen, compPosDen, intResComp, fLogDenScale defined
	dInfo
}


.checkBlock <- function(
	### filling default values in block description
	block				##<< list to be checked for entries
	,dInfos				##<< checked list of dInfo. See the argument in \code{\link{twDEMCBlockInt}}.
	,parNames			##<< parameter dimension names
){
	##details<< if fLogDenInfoId is not specified, assume that it uses the entry in fLogDenInfos	
	if( 0 == length(block$dInfoPos) ) block$dInfoPos=1
	if( 1 != length(block$dInfoPos) ) stop(".checkBlock: dInfoPos must be of length 1.")
	lpOrig <- block$dInfoPos
	if( is.character(lpOrig) ) block$dInfoPos <- match(lpOrig,names(dInfos))
	if( is.na(block$dInfoPos) ){
		stop(".checkBlock: dInfoPos specified in block are not in names of dInfos.")
	}
	dInfo <- dInfos[[block$dInfoPos]]
	
	if( 0 == length(block$compPos) ) stop(".checkBlock: entry compPos must be specified, either as index or as variable names")
	cpOrig <- block$compPos
	if( is.character(cpOrig) ){
		block$compPos <- match( cpOrig, parNames)
	}
	if( any(is.na(block$compPos)) ){
		stop(".checkBlock: some of compPos specified in block are not variable names.")
	}
	
	block$compPosInDen <- match(block$compPos,dInfo$compPosDen ) 
	if( any(is.na(block$compPosInDen)) ){
		stop(".checkBlock: some of compPos specified in block are not in compPosDen of fLogDenInfo.")
	}
	
	if( 0==length(block$fUpdateBlock) ) block$fUpdateBlock <- updateBlockTwDEMC
	if( 1!=length(block$fUpdateBlock) || !is.function(block$fUpdateBlock) ) stop(".checkBlock: fUpdateBlock must be a single function.")
	
	if( 0==length(block$requiresUpdatedDen) ) block$requiresUpdatedDen <- TRUE
	
	### block with entries compPos, dInfoPos, fUpdateBlock, requiresUpdatedDen definded
	block
}

.getResFLogDenNames <- function(
	### Extracting/Creating names for result components of logDenComp
	logDenComp		##<< vector or array (vars in rows) of results of logDensity function 
){
	##details<< 
	## Names of the result components of fLogDen are used to ditinguish internal components.
	## If the density function does not provide names, they are created.
	## If the density function has only one result component, the name is lost during parallel execution.
	## Hence, a single name is also re-created, to avoid errors on checking.
	logDenName <- "logDen"	
	if( is.array(logDenComp)){
		res <- rownames(logDenComp)
		return( if( !is.null(res) ) res else paste(logDenName,1:nrow(logDenComp),sep="") )
	}else{
		res <- names(logDenComp)
		return( if( !is.null(res) && length(res) > 1) res else paste(logDenName,1:length(logDenComp),sep="") )
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
	## \code{\link{twDEMCBlockInt}}
	## \code{\link{.generateXPropThin}}
	nChainPop = dim(Z)[3] 
	nParm = dim(Z)[2]
	##details<<  
	## Random states for chains for difference vector generation are within subsets of populations.
	## This allows simulating independent population of chains.
	## The acceptance rate may differ amonng populations. Hence, the set of previous generations to 
	## randomly select from also differs between poplations.
	# integer array (thinSteps*nChain*4) sample of chains, within populations
	X <- adrop(Z[mZ,,,drop=FALSE],1)		# assume current state X as beginning of the interval
	nStates <- nSample*nChainPop*3			# need to sample three states for snooker update
	sChains <- 1:nChainPop 
	sGens <- (mZ-nGenBack+1):mZ 
	rrGenPop <-  sample(sGens, nStates, replace = TRUE)
	rrChainsPop <-  sample(sChains, nStates, replace = TRUE)
	# in order to constrain two dimensions at the same time use the [] subset with an array see ?"["
	rr <- cbind(rrGenPop,rrChainsPop)
	#zLogDen <- array(rLogDen[rr], dim=c(1,nSample*nChainPop,3), dimnames=list(parms="logDen", steps=NULL, zi=NULL) )
	rrParms <- cbind( rep(rrGenPop,each=nParm), rep(1:nParm, nStates), rep(rrChainsPop,each=nParm) )
	zParms <- array(Z[rrParms], dim=c(nParm,nSample*nChainPop,3), dimnames=list( parms=rownames(Z), steps=NULL, zi=NULL) )
	# note that parameters are the first dimension		
	#rrParms <- cbind( rep(rrGenPop,nParm), rep(1:nParm, each=nStates), rep(rrChainsPop,nParm) )
	#zParms <- matrix( Z[rrParms], ncol=nParm)
	#z0 <- abind( zParms, zLogDen, along=1)
	
	# append as forth initial state x vector to z
	# assume that chains is the last dimension, for each step assume the same initial state x 
	chainsZ <- rep(sChains,each=nSample)
	Xs <- array(X[,chainsZ], dim=c(nParm,nSample*nChainPop,1), dimnames=list(parms=rownames(Z), steps=NULL, zi=NULL) )
	#XsLogDen <- array( rLogDen[mZ,chainsZ],dim=c(1,nSample*nChainPop,1), dimnames=list(parms="logDen", steps=NULL, zi=NULL)  ) 
	#Xs <- abind( Xs, XsLogDen, along=1)	#logDen row never used for X
	zx <- abind( zParms, Xs, along=3)
	
	##details<< 
	## States z2 is  attempted to be distinct.
	## And z3 is attempted to be distinct from x (assuming steps in z are stacked chains collectives of X)  
	iSame <- twWhichColsEqual( adrop(zx[,,2 ,drop=FALSE],3), adrop(zx[,,1 ,drop=FALSE],3) )
	i <- 0
	while( 0 < length(iSame) && i<5){
		zx[,iSame,2] <- zx[,sample.int(dim(zx)[2],length(iSame)),2]
		iSame <- twWhichColsEqual( adrop(zx[,,2 ,drop=FALSE],3), adrop(zx[,,1 ,drop=FALSE],3) )
		i<-i+1
	}
	#(tmp <- which( zx[,,2] == zx[,,1], arr.ind=TRUE))
	#when z3 is the same as x, i.e. zx[,,4]
	iSame <- twWhichColsEqual( adrop(zx[,,3 ,drop=FALSE],3), adrop(zx[,,4 ,drop=FALSE],3) )
	#iSame <- unique(which( (zx[-nrow(zx),,3]==Xs), arr.ind=TRUE )[,2])
	i <- 0
	while( 0 < length(iSame) && i<5){
		zx[,iSame,3] <- zx[,sample.int(dim(zx)[2],length(iSame)),3]
		iSame <- twWhichColsEqual( adrop(zx[,,3 ,drop=FALSE],3), adrop(zx[,,4 ,drop=FALSE],3) )
		i<-i+1
	}
	#(tmp <- which( zx[,,3] == zx[,,4], arr.ind=TRUE))
	### random states (Nparms+1,(steps*nChainPop), 4)
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
	## \code{\link{twDEMCBlockInt}}
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
	xStepAndExtra <- resArraySteps <- array(res, dim=c(nParm+1,nSample,nChain), dimnames=list(parms=rownames(res),steps=NULL,chains=NULL) )
	#deprecated expected steps as last dimension
	#xStepAndExtra <- aperm(resArraySteps,c(1,3,2))
	# numeric array (Npar+1,Nchains,Nsteps): difference vectors in parameter space for steps and chains
	# last row is the extra LogDen associated with snooker update
	
	xStep <- xStepAndExtra[-nrow(xStepAndExtra),,,drop=FALSE]			#dito
	rExtra <- adrop(xStepAndExtra[nrow(xStepAndExtra),,,drop=FALSE],1)			#second dim (columns step within Thinning interval)
	list(xStep=xStep, rExtra=rExtra)
	### List with components \describe{
	### \item{xStep}{numeric array (Npar x NSteps x Nchain): difference vectors in parameter space}
	### \item{rExtra}{numeric matrix (Nsteps x NChain): some extra LogDensity from snooker update}}
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

.updateIntervalTwDEMCParChains <- function(
        ### Perform update steps within next thinning interval parellelizing by chains.
        # first three arguments relate to current state (stored in chainState in twDEMCBlockInt)
        X				##<< numeric matrix: current location of chains rows: parameters, columns: chains
        ,logDenCompX	##<< numeric matrix: result of calls to density functions from previous step
        ,parUpdateDen	##<< for each parameter/density combination: is the density up to date from previous step
        ,xStep 			##<< array rows: difference vectors in parameter space, cols: steps, zdim: chains
        ,rExtra			##<< numeric matrix: row: step within thinning interval, col: chain
        ,tempDenCompSteps	##<< numeric array:  (nGenThin x nResComp x nsPop)
        ,nsPop			##<< the number of populations that take part in this step
        ,pAcceptPar		##<< numeric matrix: current acceptance rate for each (block x population)
        #argsUpdateBlocksTwDEMC, ##<< list of arguments to updateBlocksTwDEMC
        #argsDEMCStep,	##<< see \code{\link{.doDEMCStep}}
        ,remoteFunArgs	
        ### see \code{\link[twSnowfall]{sfRemoteWrapper}}
        ### Must include entries \itemize{
        ### \item{remoteFun=.updateBlocksTwDEMC}
        ### \item{argsUpdateBlocksTwDEMC (maybe character name of exported variable}
        ### \item{remoteDumpfileBasename} 
        ### }
        ### Can be a name of a previously exported variable.
        ,debugSequential=FALSE	##<< see \code{\link[twSnowfall]{sfFArgsApplyDep}}
        ,freeMasterNode=FALSE	##<< if set to TRUE, no job is submitted to first node, so that this node can dispatch jobs without waiting, see \code{\link[twSnowfall]{sfFArgsApplyDep}}
        ,isRecordProposalsPop=rep(FALSE,nsPop)	##<< logical vector: for each population: if TRUE then proposals and results of rLogDen are recorded in result$Y.
        , chainsPop=matrix(1:ncol(X),ncol=nsPop, dimnames=list(iChainInPop=NULL, iPop=NULL))		##<< list of integer vectors:
        , isPops			##<< integer vector of populations that have not yet dropped out
){
    # .updateIntervalTwDEMCParChains
    dimXStep <- dim(xStep)
    nStep <- dimXStep[2]
    nsChain <- dimXStep[3]
    nChainPop = nsChain / nsPop
    #
    #iChain <- 1
    fapply <- if( isTRUE(debugSequential) ) lapply else sfLapply
    res <- fapply( 1:nsChain, function(iChain){
                iPop <- ((iChain-1) %/% nChainPop) +1  
                .updateIntervalTwDEMCChain( X[,iChain], logDenCompX[,iChain]
                    , adrop(parUpdateDen[,,iChain,drop=FALSE],3)
                    , adrop(xStep[,,iChain,drop=FALSE],3)
                    , rExtra[,iChain]
                    , adrop(tempDenCompSteps[,,iPop,drop=FALSE],3)
                    , pAcceptPar[,iPop]
                    ,remoteFunArgs=remoteFunArgs
                    ,isRecordProposal[iPop]
                    , iPop = iPop
                )                
            }) 
    ##value<< list with components
    res1 <- res[[1]]
    X <- array( NA_real_,  c(length(res1$X), nsChain), dimnames=list(parms=names(res1$X), iChain=NULL))
    logDenCompX <- array( NA_real_,  c(length(res1$logDenCompX), nsChain), dimnames=list(parms=names(res1$logDenCompX), iChain=NULL))
    parUpdateDen <- array( NA,  c(dim(res1$parUpdateDen), nsChain), dimnames=c(dimnames(res1$parUpdateDen), iChain=NULL))
    accepted <- array( NA_real_,  c(length(res1$accepted), nsChain), dimnames=list(iBlock=names(res1$accepted), iChain=NULL))
    Y <-  array( NA_real_,  c(dim(res1$Y), nsChain), dimnames=c(dimnames(res1$Y), iChain=NULL))
    for( iChain in 1:nsChain ){
        resi <- res[[iChain]]
        X[,iChain] <- resi$X
        logDenCompX[,iChain] <- resi$logDenCompX
        parUpdateDen[,,iChain] <- resi$parUpdateDen
        accepted[,iChain] <- resi$accepted
        Y[,,iChain] <- resi$Y
    } 
    resDo <- list(	##describe<< 
            X=X							##<< matrix current position, column for each chain
            , logDenCompX=logDenCompX	##<< matrix: result components of fLogDen current position, column for each chain 
            #, logDenX=logDenX			##<< vector current logDen of chains
            , parUpdateDen=parUpdateDen	##<< numeric matrix (par x den x chain): for each parameter/density combination: is the density up to date
            , accepted=accepted			##<< numeric matrix (blocks x chains): number of accepted steps for each chain
            , Y=if( !is.null(res1$Y) ) Y else NULL	##<< numeric matrix (steps x 2+params+result x chains): accepted, rLogDen, parms, and all fLogDen result components for each proposal
    ) ##end<< 
}

.updateIntervalTwDEMCChain <- function(
        ### Perform update steps within next thinning interval for a single chain (called remotely from .updateIntervalTwDEMCParChains) 
        # first three arguments relate to current state (stored in chainState in twDEMCBlockInt)
        X				##<< numeric vector: current location
        ,logDenCompX	##<< numeric vector: result of calls to density functions from previous step
        ,parUpdateDen	##<< for each parameter/density combination: is the density up to date from previous step
        ,xStep 			##<< matrix rows: difference vectors in parameter space, cols: steps
        ,rExtra			##<< numeric vector: row: step within thinning interval
        ,tempDenCompSteps	##<< numeric matrix:  (nGenThin x nResComp)
        #,nsPop			##<< the number of populations that take part in this step
        ,pAcceptPar		##<< numeric vector: current acceptance rate for each (block)
        #argsUpdateBlocksTwDEMC, ##<< list of arguments to updateBlocksTwDEMC
        #argsDEMCStep,	##<< see \code{\link{.doDEMCStep}}
        ,remoteFunArgs	
        ### see \code{\link[twSnowfall]{sfRemoteWrapper}}
        ### Must include entries \itemize{
        ### \item{remoteFun=.updateBlocksTwDEMC}
        ### \item{argsUpdateBlocksTwDEMC (maybe character name of exported variable}
        ### \item{remoteDumpfileBasename} 
        ### }
        ### Can be a name of a previously exported variable.
        ,debugSequential=FALSE	##<< see \code{\link[twSnowfall]{sfFArgsApplyDep}}
        #,freeMasterNode=FALSE	##<< if set to TRUE, no job is submitted to first node, so that this node can dispatch jobs without waiting, see \code{\link[twSnowfall]{sfFArgsApplyDep}}
        ,isRecordProposalsPop=FALSE 	##<< logical
        , iPop          ##<< the population that the current chain is in
){
    # .updateIntervalTwDEMCPar
    ##seealso<< 
    ## \code{\link{twDEMCBlockInt}}
    ## \code{\link{.doDEMCStep}}
    dimXStep <- dim(xStep)
    nStep <- dimXStep[2]
    #d <- as.list(structure(dimXStep,names=c("parms","chains","steps"))) 
    iGenT <- (1:nStep)
    nCases = nStep
    if( !(nCases == length(rExtra)) )
        stop("number of cases in steps must correspond to length of rExtra")
    res0 <- list(
                        xC=X
                        ,logDenCompC=logDenCompX
                        ,parUpdateDenC=parUpdateDen
                )
    #
    # record the proposals, their density results and their acceptance for each generation in thinning interval
    parNames <- names(X)
    resCompNames <- names(logDenCompX)
    if( is.null(resCompNames) ) resCompNames <- paste("lDen",seq_along(res[[1]]$logDenCompP),sep="") 
    tmp <- c(parNames,paste("accepted",seq_along(pAcceptPar),sep=""),resCompNames)
    Y <- array( NA_real_, dim=c(nStep=nStep,comp=length(tmp)), dimnames=list(steps=NULL,comp=tmp) )			
    #
    # recored acceptance state
    acceptedSteps <- array( NA_real_, dim=c(nStep=nStep,block=length(pAcceptPar)) )
    #
    #i <- 1
    prevRes <- res0
    tmpf <- remoteFunArgs$argsUpdateBlocksTwDEMC
    for( iStep in  1:nStep){
        args <-c( prevRes[1:3] 
                ,list( step=xStep[,iStep] 
                        , rExtra=rExtra[iStep ,drop=TRUE]
                        #, tempGlobalC=tempGlobalSteps[iStep0+1,isPop ,drop=TRUE]
                        , tempDenCompC=tempDenCompSteps[iStep,, drop=TRUE]
                        , pAccept=pAcceptPar
                        #, isPop=isPop
                        , iPop=iPop
                    )
                    , within(remoteFunArgs, remoteFun <- NULL)
                    )
         prevRes <- do.call( .updateBlocksTwDEMC, args )
         Y[iStep,] <- c( prevRes$xP, prevRes$accepted, prevRes$logDenCompP )
         acceptedSteps[iStep,] <- prevRes$accepted
         prevRes$accepted
     }
     #
    ##value<< list with components
    resDo <- list(	##describe<< 
            X=prevRes$xC			    ##<< vector (nParm) current position
            , logDenCompX=prevRes$logDenCompC	##<< vector (nResComp): result components of fLogDen current position 
            #, logDenX=logDenX			##<< vector current logDen of chains
            , parUpdateDen=prevRes$parUpdateDenC	##<< numeric matrix (par x den): for each parameter/density combination: is the density up to date
            , accepted=colSums(acceptedSteps)			##<< numeric vector (blocks): number of accepted steps for each chain
            , Y=Y 						##<< numeric matrix (steps x 2+params+result): accepted, rLogDen, parms, and all fLogDen result components for each proposal
    ) ##end<< 
}


.updateIntervalTwDEMCPar <- function(
	### Perform update steps within next thinning interval for all the chains in load balanced parallel way.
	# first three arguments relate to current state (stored in chainState in twDEMCBlockInt)
	X				##<< numeric matrix: current location of chains rows: parameters, columns: chains
	,logDenCompX	##<< numeric matrix: result of calls to density functions from previous step
	,parUpdateDen	##<< for each parameter/density combination: is the density up to date from previous step
	,xStep 			##<< array rows: difference vectors in parameter space, cols: steps, zdim: chains
	,rExtra			##<< numeric matrix: row: step within thinning interval, col: chain
	,tempDenCompSteps	##<< numeric array:  (nGenThin x nResComp x nsPop)
	,nsPop			##<< the number of populations that take part in this step
	,pAcceptPar		##<< numeric matrix: current acceptance rate for each (block x population)
	#argsUpdateBlocksTwDEMC, ##<< list of arguments to updateBlocksTwDEMC
	#argsDEMCStep,	##<< see \code{\link{.doDEMCStep}}
	,remoteFunArgs	
	### see \code{\link[twSnowfall]{sfRemoteWrapper}}
	### Must include entries \itemize{
	### \item{remoteFun=.updateBlocksTwDEMC}
	### \item{argsUpdateBlocksTwDEMC (maybe character name of exported variable}
	### \item{remoteDumpfileBasename} 
	### }
	### Can be a name of a previously exported variable.
	,debugSequential=FALSE	##<< see \code{\link[twSnowfall]{sfFArgsApplyDep}}
	,freeMasterNode=FALSE	##<< if set to TRUE, no job is submitted to first node, so that this node can dispatch jobs without waiting, see \code{\link[twSnowfall]{sfFArgsApplyDep}}
	,isRecordProposalsPop=rep(FALSE,nsPop)	##<< logical vector: for each population: if TRUE then proposals and results of rLogDen are recorded in result$Y.
	, chainsPop=matrix(1:ncol(X),ncol=nsPop, dimnames=list(iChainInPop=NULL, iPop=NULL))		##<< list of integer vectors:
	, isPops			##<< integer vector of populations that have not yet dropped out
){
	# .updateIntervalTwDEMCPar
	##seealso<< 
	## \code{\link{twDEMCBlockInt}}
	## \code{\link{.doDEMCStep}}
	
	##details<< 
	## The step must be the last dimension in all arguments in order to make use of dependence step 
	## in load balanced execution.
	
	#if(!is.numeric(X) | !is.numeric(logDenCompX) | !is.numeric(logDenX) )
	#	stop(".doDEMCSteps: first three arguments must be numeric")
	#.tmpf <- function(iStep){adrop(xStep[,iStep, ,drop=FALSE],2)}	# all chains of current step 
	#xStepStacked <- do.call( rbind, lapply(1:(dim(xStep)[3]),.tmpf) )
	#xStepStacked <- abind( lapply(1:(dim(xStep)[2]),.tmpf), along=2 )
	#all chains of first step, all chains of second step, ...
	
	dimXStep <- dim(xStep)
	nsChain <- dimXStep[3]
	nStep <- dimXStep[2]
	#d <- as.list(structure(dimXStep,names=c("parms","chains","steps"))) 
	iGenT <- (1:nStep)
	nCases = nsChain * nStep
	nChainPop = nsChain / nsPop
	if( !(nCases == length(rExtra)) )
		stop("number of cases in steps must correspond to length of rExtra")
	#iPops <- matrix( 1:nChain, ncol=nPops)	
	#F_APPLY <- .updateBlocksTwDEMC	#the function executed on the nodes: one cycle of block updates
	F_APPLY <- sfRemoteWrapper		# reuse the exported arguments in remoteFunArgs, which is a name	
	#fArgsSame <- list( remoteFunArgs=as.name("argsDEMCStepWrapper") )	#exported, includes remoteFun=.doDEMCStep, remoteDumpfileBasename and argsDEMCStep
	F_ARGS <- function(i,prevRes){ 
		iStep0 <- ((i-1) %/% nsChain) # i is the state processed, c( iChains_1Step, iChains_2Step ... )
		iChain0 <- (i-1) - iStep0*nsChain   #((i-1) %% nChain)	 
		isPop <-(iChain0 %/% nChainPop)+1
		#args <-c( prevRes[c("xC", "logDenCompC", "parUpdateDenC")],
		args <-c( prevRes[1:3] 
			,list( step=xStep[,1+iStep0,1+iChain0 ,drop=TRUE]
				, rExtra=rExtra[1+iStep0,1+iChain0 ,drop=TRUE]
				#, tempGlobalC=tempGlobalSteps[iStep0+1,isPop ,drop=TRUE]
				, tempDenCompC=tempDenCompSteps[iStep0+1,,isPop ,drop=TRUE]
				, pAccept=pAcceptPar[,isPop]
				#, isPop=isPop
				, iPop=isPops[isPop]
				, remoteFunArgs=remoteFunArgs )
		)}	
	#mtrace(F_ARGS)
	#.res0 <- lapply(1:nrow(X),function(row){X[row,]})
	res0 <- lapply(1:nsChain,function(iChain){list(
				xC=X[,iChain,drop=TRUE]
				,logDenCompC=logDenCompX[,iChain,drop=TRUE]
				,parUpdateDenC=adrop(parUpdateDen[,,iChain,drop=FALSE],3)
			)})
	# all chains of first step, all chains of second step, ...
	res <- sfFArgsApplyDep( nCases, F_ARGS, F_APPLY, res0, debugSequential=debugSequential)
	# extract the state of the last step
	# modify in place, so that dimnames etc are preserved
	endChain0 <- (nStep-1)*nsChain	#index before last step of all chains
	for( iChain in 1:nsChain ){
		resChain <- res[[endChain0+iChain]]
		X[,iChain] <- resChain$xC
		logDenCompX[,iChain] <- resChain$logDenCompC
		parUpdateDen[,,iChain] <- resChain$parUpdateDenC
	}
	
	# record the proposals, their density results and their acceptance for each generation in thinning interval
	Y <- NULL
	if( any(isRecordProposalsPop) ){
		parNames <- rownames(X)
		resCompNames <- names(res[[1]]$logDenCompP)
		if( is.null(resCompNames) ) resCompNames <- paste("lDen",seq_along(res[[1]]$logDenCompP),sep="") 
		tmp <- c(parNames,paste("accepted",seq_along(res[[1]]$accepted),sep=""),resCompNames)
		Y <- array( NA_real_, dim=c(nStep=nStep,comp=length(tmp),nChain=nsChain)			, dimnames=list(steps=NULL,comp=tmp,chains=NULL) )			
		for( iPop in which(isRecordProposalsPop)){
			for(iStep in 1:nStep){
				chain0 <- (iStep-1) * nsChain	# position before first chain of step iStep
				for( iChain in chainsPop[,iPop]){
					i <- chain0 +iChain
					Y[iStep,,iChain] <- c(	res[[i]]$xP, res[[i]]$accepted, res[[i]]$logDenCompP )
				} # iChain
			} # iStep
		} #iPop
	}
	
	# record acceptance rate
	# need cbind because sapply will produce a vector instead of a matrix for only one block
	acceptedStates <- do.call( cbind, lapply( res, "[[", "accepted" )) # need to be stacked by chains but steps is last dim 
	acceptedM <- array(acceptedStates, dim=c(nrow(acceptedStates),nsChain,nStep) ) # blocks x  chains x steps 	
	accepted <- apply(acceptedM,c(1,2),sum)
	dimnames(accepted) = list(blocks=NULL, chains=NULL)
	#accepted <- matrix(0.25*remoteFunArgs$argsUpdateBlocksTwDEMC$argsFUpdateBlocks[[1]]$ctrl$thin	, nrow=length(res[[1]]$accepted), ncol=nChain )
	
	##value<< list with components
	resDo <- list(	##describe<< 
		X=X							##<< matrix current position, column for each chain
		, logDenCompX=logDenCompX	##<< matrix: result components of fLogDen current position, column for each chain 
		#, logDenX=logDenX			##<< vector current logDen of chains
		, parUpdateDen=parUpdateDen	##<< numeric matrix (par x den x chain): for each parameter/density combination: is the density up to date
		, accepted=accepted			##<< numeric matrix (blocks x chains): number of accepted steps for each chain
		, Y=Y 						##<< numeric matrix (steps x 2+params+result): accepted, rLogDen, parms, and all fLogDen result components for each proposal
	) ##end<< 
}

.updateBlocksTwDEMC <- function(
	### One step of updating all blocks for one chain (called remotely)
	xC				##<< current state
	,logDenCompC	##<< result of calls to density functions from previous step
	,parUpdateDenC	##<< for each parameter/density combination: is the density up to date from previous step
	,step			##<< proposed step
	,rExtra			##<< adjustment of Metropolis state
	#,tempGlobalC			##<< current global temperature for given population
	,tempDenCompC	##<< numeric vector of length(logDenCompC): current temperature for each result component
	,pAcceptC		##<< numeric vector (nBlock) current acceptance rates of the blocks
	,iPop			##<< population in total array of populations
	,argsUpdateBlocksTwDEMC	##<< further arguments that do not change between calls (exported and handled by remoteWrapper)
){
	if( !all(is.finite(pAcceptC)))
		stop(".updateBlocksTwDEMC: encountered non-finited pAcceptC")
	a <- argsUpdateBlocksTwDEMC
	#with( argsUpdateBlocksTwDEMC,{
	acceptedBlock <- numeric(a$nBlock) 
	xP <- xC		# place for recording proposal
	logDenCompP <- logDenCompC # place for recording density results for proposal
    #
    intermediates <- list()     # initially empty
	#
    #recover()
	#iBlock<-1
    #iBlock<-2
	#dInfo1 <- a$argsFUpdateBlocks[[1]]
	for( iBlock in seq_along(a$argsFUpdateBlocks)){
		#block <- blocks[[iBlock]]
		blockArgs <- a$argsFUpdateBlocks[[iBlock]] #merged the information in argsFUpdateBlocks
        usesIntermediate <- length(blockArgs$intermediate) != 0
        #
		cd <- blockArgs$compPosDen # components used by the current density function
		cu <- blockArgs$compPos	   # components to be updated in this block
		#cuInD <- block$compPosInD 
		##details<< \describe{}{\item{Recalculating logDensity of accepted state.}{
		## If one of the components used by the current density function has been
		## updated against another density function, then the current density of  
		## the accepted state is out of date, and needs to be recalculated.
		##}}
		if( blockArgs$requiresUpdatedDen && !all( parUpdateDenC[blockArgs$dInfoPos,cd] ) ){
            # intermediate state may have been calculated
            blockArgs$argsFLogDen$intermediate <- if( usesIntermediate )  intermediates[[blockArgs$intermediate]] else NULL
            # here do not allow for internal rejection
            logDenCompC[blockArgs$resCompPos] <- tmp <- blockArgs$fLogDenScale * do.call( blockArgs$fLogDen, c(list(xC[cd]), blockArgs$argsFLogDen) )	# evaluate logDen
			parUpdateDenC[blockArgs$dInfoPos,cd] <- TRUE
            if( usesIntermediate ){
                intermediates[[ blockArgs$intermediate]] <- attr(tmp,"intermediate") 
            } 
		}	
		# make sure to use names different from already existing in argsFUpdateBlocks[[iBlock]]
		# in c latter entries overwrite previous entries, so use argsFUpdateBlocks first
        blockArgs$argsFLogDen$intermediate <- NULL      # itermediate refers to current state not for proposal
		argsFUpdateBlock <- c( blockArgs, list(
				step=step
				, rExtra=rExtra
				, logDenCompC=logDenCompC[ blockArgs$resCompPos  ]
				#, tempC=tempGlobalC
				, tempDenCompC=tempDenCompC[ blockArgs$resCompPos ]	
				, pAccept=pAcceptC
				, upperParBounds=a$upperParBoundsPop[[iPop]]
				, lowerParBounds=a$lowerParBoundsPop[[iPop]]
				, iPop=iPop
			))
		#resUpdate <- updateBlockTwDEMC( xC[cd], argsFUpdateBlock=argsFUpdateBlock )
		resUpdate <- blockArgs$fUpdateBlock( xC[cd], argsFUpdateBlock=argsFUpdateBlock )       # call the update
		acceptedBlock[iBlock] <- resUpdate$accepted # record acceptance
		if( !all(is.finite(resUpdate$accepted)))
			stop(".updateBlocksTwDEMC: encountered non-finited acceptance rate")
		if( resUpdate$accepted != 0){
			#if( resUpdate$xC["a"] > 10.8) recover()
			xC[cu] <- xP[cu] <- resUpdate$xC		# update the components of current state
            #undebug(setIntermediate)
			parUpdateDenC[,cu] <- FALSE	# indicate all density results out of date for the update parameters
			if( 0 != length(resUpdate$logDenCompC) ){				
                if( usesIntermediate) 
                    intermediates[[ blockArgs$intermediate]] <- attr(resUpdate$logDenCompC,"intermediate")  # update the intermediate
                logDenCompC[ blockArgs$resCompPos  ] <- logDenCompP[ blockArgs$resCompPos  ] <- resUpdate$logDenCompC # update the result of the used density function
				parUpdateDenC[blockArgs$dInfoPos,cu] <- TRUE	# if density has been calculated, indicate that it is up to date 
            }
		}else{ # not accepted
			if( 0 != length(resUpdate$xP) ) xP[ cu  ] <- resUpdate$xP 				# proposal 		
			if( 0 != length(resUpdate$logDenCompP) ) logDenCompP[ blockArgs$resCompPos  ] <- resUpdate$logDenCompP # density results for proposal
            # keep with previous intermediate                
		}
	} #iBlock

	##value<< list with components
	resUpdate <- list( ##describe <<
		xC=xC						##<< numeric vecotr: current accepted state
		,logDenCompC=logDenCompC	##<< numeric vector: result components of fLogDen for current position
		,parUpdateDenC=parUpdateDenC			##<< integer vector: logDensity that recently updated parameter at given index
		,accepted=acceptedBlock		##<< numeric vector: acceptance of each block (0: not, 1: did, (0..1): DRstep)
		,xP=xP						##<< numeric vector: result components of fLogDen for proposal
        #,intermediate
		,logDenCompP=logDenCompP
	) ##enf<<
	#})	# with
}

### \item{tempC}{global temperature}
### \item{TFix}{named numeric vector (comp): components with fixed Temperature}
### \item{posTFix}{integer vector (comp): =match(TFix, compNames): positions of TFix within comp provided for performance reasons}
### \item{ctrl$useMultiT}{boolean wheter to scale temperature for different data streams}
### \item{TProp}{numeric vector of length(resCompPos): proportions of temperature for different data streams}

#undebug(updateBlockTwDEMC)
updateBlockTwDEMC <- function( 
	### Perfrom one DEMC step, function to be called in remote process.
	xC				##<< numeric vector: current state that is used in density function of the block
	,argsFUpdateBlock	
### further argument provided for generating the update  \describe{
### \item{compPosInDen}{ positions of the dimensions in x that are updated in this block}
### \item{step}{proposed jump}
### \item{rExtra}{extra portion in Metropolis decision because of selecting the jump}
### \item{logDenCompC}{numeric vector: former result of call to same fLogDen}
### \item{tempDenCompC}{numeric vector of length(logDenCompC): temperature for each density result component}
### \item{fDiscrProp,argsFDiscrProp}{function and additional arguments applied to xProp, e.g. to round it to discrete values}
### \item{argsFLogDen, fLogDenScale}{additional arguments to fLogDen and scalar factor applied to result of fLogDen}
### \item{posLogDenInt}{the matching positions of intResCompNames within the the results components that are handled internally}
### \item{ctrl$DRgamma}{ if !0 and >0 delayed Rejection (DR) (Haario06) is applied by jumping only DRgamma distance along the proposal }
### \item{upperParBounds}{ named numeric vector, see \code{\link{twDEMCBlockInt}}  }
### \item{lowerParBounds}{ named numeric vector, see \code{\link{twDEMCBlockInt}}  }
### \item{fCalcComponentTemp}{ functiont to calculate temperature of result components, (way of transporting calcComponentTemp to remote process) }
### \item{iPop}{ just for debugging }
### }
){
    #recover()
	#updateBlockTwDEMC
	##seealso<< 
	## \code{\link{twDEMCBlockInt}}
	## \code{.updateBlocksTwDEMC}
	# # \code{\link{.updateBlocksTwDEMC}}
	#attach( argsDEMCStep )
	#stop(".doDEMCStep: stop to trace error in remote R invocation.")
	a <- argsFUpdateBlock
	#with( argsFUpdateBlock, {		# with produces significant performance loss
	boResFLogDenX <- (length(a$intResCompPos) > 0)
	# LogDensity of accepted state
	#La <- logDenCompC	#logDensity	components of accepted state
	if( 0 == length(a$logDenCompC)) stop("updateBlockTwDEMC: encountered empty logDenCompC")			
	#assume that all is.finite(logDenCompC), make sure in twDEMCBlockInt
	LaExt <- La <- a$logDenCompC
	#logDenC <- sum(La)
	#a$tempDenCompC <- a$fCalcComponentTemp( temp=a$tempC, TFix=a$TFix, TProp=a$TProp, useMultiT=a$ctrl$useMultiT,posTFix=a$posTFix)
	#
	accepted<-FALSE
	xP <- xC
	xP[a$compPosInDen] <- xP[a$compPosInDen] + a$step[a$compPosInDen]  
	#
	boOutside <- 
		any( sapply( names(a$upperParBounds), function(pname){ xP[pname] >= a$upperParBounds[pname] })) ||
		any( sapply( names(a$lowerParBounds), function(pname){ xP[pname] < a$lowerParBounds[pname] }))
	#if(xProp[1] < 10.8)	recover()
	#sapply( names(a$lowerParBounds), function(pname){ pname })
	#sapply( names(a$lowerParBounds), function(pname){ xP[pname] })
	#sapply( names(a$lowerParBounds), function(pname){ a$lowerParBounds[pname] })
	if( boOutside ){
		# bp(11)
		# if it is still outside (maybe opposite border) reject step and give -Inf as logDenResult
		#logDenP=logAlpha10=-Inf		
		logAlpha10=-Inf		#logAlpha10 is log of the initial acceptance ratio for DR step (-Inf no chance of acceptance)
		Lp <- a$logDenCompC	#results for the proposal
		Lp[] <- -Inf
	}else{
		# discrtize proposal
		if( is.function(a$fDiscrProp)) xP = do.call(a$fDiscrProp,xP,a$argsFDiscrProp, quote=TRUE)
		Lp <- a$fLogDenScale * if(boResFLogDenX){
				a$logDenCompIntWithin <- a$logDenCompC[a$intResCompPos]
				TiInt <- a$tempDenCompC[a$intResCompPos]
				do.call( a$fLogDen, c(list(xP, a$logDenCompIntWithin, TiInt), a$argsFLogDen) )	# evaluate logDen
			}else{
				do.call( a$fLogDen, c(list(xP), a$argsFLogDen) )	# evaluate logDen
				#a$fLogDen( xProp, blockIndices=a$argsFLogDen$blockIndices, fModel=a$argsFLogDen$fModel, obs=a$argsFLogDen$obs, invCovar=a$argsFLogDen$invCovar, xval=a$argsFLogDen$xval )
			}
		#take care that the result has always the same sames, even when if fails
		#if( 0==length(names(res)))
		#	stop("encountered result of fLogDen without names")
		#if( !identical(names(logDenCompC),names(res)))
		#	stop("encountered result with different names")
		#strip attributes other than names, else twDynamicClusterApplyDep fails with big data chunks
		#attributes(Lp) <- list(names=names(Lp))
		#logDenP=-Inf
		#make sure Lp, La have the same order and legnth
		#if( !identical( names(Lp), names(La)) ) stop(".doDEMCStep: logDenCompC must contain the same components and the order of result of fLogDen." )
		if( all(is.finite(Lp))){
			#logDenP <- sum(Lp)
			##details<< \describe{\item{internal Metropolis step}{
			## if posLogDenInt is given, then these components of result of fLogDen are handled
			## internally. Hence, for Metropolis step here operates only on other remaining components.
			## }}
			posTExt <- setdiff( seq_along(a$resCompPos), a$intResCompPos )		#externally handled components
			nExt <- length(posTExt)
			#posTFixExt <- setdiff(posTFix,posLogDenInt)		#externally handled components with fixed temperature
			#posTVarExt <- setdiff(seq_along(Lp), c(posTFix,posLogDenInt))	#externally handled componetns with variable temperature
			#nFixExt <- length(posTFixExt)
			#nVarExt <- length(posTVarExt)
			#nExt <- nFixExt + nVarExt
			
			#Metropolis step
			#logr = (logDenPropExt+rExtra - logDenXExt) / tempC
			logrDS10 <- (Lp[posTExt]-La[posTExt])/a$tempDenCompC[posTExt]
			logAlpha10 <- a$rExtra + sum(logrDS10) 
			accepted <-  if( is.finite(logAlpha10) && (logAlpha10  > log(runif(1)) ) ) 1 else 0
			if(accepted != 0){
				xC <- xP
				#logDenC <- logDenP
				a$logDenCompC <- Lp
			}				
		}else logAlpha10 <- -Inf
	} # end check outside parBounds
	
	if(!accepted && !is.null(a$ctrl$DRgamma) && (a$ctrl$DRgamma > 0) && 
		( boOutside ||	(!is.null(a$ctrl$minPCompAcceptTempDecr) && (a$pAccept < 1.2*a$ctrl$minPCompAcceptTempDecr)))
		) {
		#----- delayed rejection (DR) step
		# only if across parBoundEdge or acceptance rate drops below 1.2*minAcceptrate
		# repeat all above with delayed rejection (DR) step, only adjust DRfac after calculating Lp
		Lp1 <- Lp
		xProp1 <- xP	# former proposal
		xP <- xC
		xP[a$compPosInDen] <- xP[a$compPosInDen] + a$ctrl$DRgamma*a$step[a$compPosInDen]				
		
		boOutside <- 
			any( sapply( names(a$upperParBounds), function(pname){ xP[pname] > a$upperParBounds[pname] })) ||
			any( sapply( names(a$lowerParBounds), function(pname){ xP[pname] < a$lowerParBounds[pname] }))
		if( !boOutside ){
			if( is.function(a$fDiscrProp)) xP = do.call(a$fDiscrProp,xP,a$argsFDiscrProp, quote=TRUE)
			Lp <- a$fLogDenScale * if(boResFLogDenX){
					a$logDenCompIntWithin <- a$logDenCompC[a$intResCompPos]
					TiInt <- a$tempDenCompC[a$intResCompPos]
					do.call( a$fLogDen, c(list(xP, a$logDenCompIntWithin, TiInt), a$argsFLogDen) )	# evaluate logDen
				}else{
					do.call( a$fLogDen, c(list(xP), a$argsFLogDen) )	# evaluate logDen
					#a$fLogDen( xProp, blockIndices=a$argsFLogDen$blockIndices, fModel=a$argsFLogDen$fModel, obs=a$argsFLogDen$obs, invCovar=a$argsFLogDen$invCovar, xval=a$argsFLogDen$xval )
				}
			#take care that the result has always the same sames, even when if fails
			#if( 0==length(names(res)))
			#	stop("encountered result of fLogDen without names")
			#if( !identical(names(logDenCompC),names(res)))
			#	stop("encountered result with different names")
			#strip attributes other than names, else twDynamicClusterApplyDep fails with big data chunks
			attributes(Lp) <- list(names=names(Lp))
			#logDenP=-Inf
			#make sure Lp, La have the same order and legnth
			#if( !identical( names(Lp), names(La)) ) stop(".doDEMCStep: logDenCompC must contain the same components and the order of result of fLogDen." )
			if( all(is.finite(Lp))){
				#logDenP <- sum(Lp)
				##details<< \describe{\item{internal Metropolis step}{
				## if posLogDenInt is given, then these components of result of fLogDen are handled
				## internally. Hence, for Metropolis step here operates only on other remaining components.
				## }}
				posTExt <- setdiff( seq_along(a$tempDenCompC), a$intResCompPos )		#externally handled components
				nExt <- length(posTExt)
				#posTFixExt <- setdiff(posTFix,posLogDenInt)		#externally handled components with fixed temperature
				#posTVarExt <- setdiff(seq_along(Lp), c(posTFix,posLogDenInt))	#externally handled componetns with variable temperature
				#nFixExt <- length(posTFixExt)
				#nVarExt <- length(posTVarExt)
				#nExt <- nFixExt + nVarExt
				#Metropolis step 
				#logr = (logDenPropExt+rExtra - logDenXExt) / tempC
				logrDS20 <- (Lp[posTExt]-La[posTExt])/a$tempDenCompC[posTExt]
				logAlpha20 <- a$rExtra + sum(logrDS20)
				#---  here correct with first stage DR factor (1-alpha21)/(1-alpha10) with meaning 0:accepted 1:first proposal 2:second proposal
				logrDS21 <- (Lp1[posTExt]-Lp[posTExt])/a$tempDenCompC[posTExt]
				logAlpha21 <- sum(logrDS21)
				logAlpha2 <- suppressWarnings( logAlpha20  +log(1-exp(logAlpha21)) -log(1-exp(logAlpha10)) )	# log and exp may produce NaNs 
				accepted <-  if( is.finite(logAlpha2) && ( logAlpha2 > log(runif(1)) ) ) a$ctrl$DRgamma else 0  
				if(accepted != 0){
					xC <- xP
					#logDenC <- logDenP
					a$logDenCompC <- Lp
				}				
			}
		} # end !boOutside in DR step 
	}	# end DR step
	
	#will invoke prevRes[c("x", "logDenCompC", "logDenAcc")]
	if(!is.numeric(xC) | !is.numeric(a$logDenCompC) ) #| !is.numeric(logDenC) )
		stop(".twDEMCBlocks: x, logDenAcc and logDenCompC must be numeric")
	
	##value<< list with components
	list(	##describe<<
		accepted=accepted			##<< boolean scalar: if step was accepted
		, xC=xC[a$compPosInDen]		##<< numeric vector: components of position in parameter space that are being updated (if accepted then the same as xP)
		, logDenCompC=a$logDenCompC	##<< numeric vector: result components of fLogDen for current position (if accepted then the same as logDenCompP)
		#, logDenC =logDenC			##<< numeric vector: summed fLogDen for current accepted position
		, xP=xP[a$compPosInDen]		##<< numeric vector: components of proposal that are being 
		, logDenCompP=Lp			##<< numeric vector: result components of fLogDen for proposal, may inlcude attr "intermediate"
	#, logDenP=logDenP			##<< numeric vector: summed fLogDen for proposal
	)	##end<<
	#})
}

setMethodS3("twDEMCBlock","array", function( 
		### Initialize \code{\link{twDEMCBlockInt}} by array of initial population and remove those generations from results afterwards
		x ##<< initial population: a numeric array (M0 x d x nChain) see details in \code{\link{twDEMCBlockInt}}  
		,...	##<< further arguments to \code{\link{twDEMCBlockInt}}
		,nPop = 1	##<< number of populations in x
		#,TProp = NULL	##<< numeric vector (nResComp): proportions of temperature for result components: equal for all populations
		, X=NULL		##<< initial state (nParm x nChain)
		, logDenCompX = NULL ##<< numeric matrix (nResComp x nChains) initial state
		, upperParBounds = vector("list",nPop)
		### list of named numeric vectors: giving upper parameter bounds for each population 
		### for exploring subspaces of the limiting distribution, see details
		### \cr Alternatively a single numeric vector can be supplied, which is replicated for each population.
		, lowerParBounds = vector("list",nPop)
	### similar to upperParBounds
	){
		#twDEMC.array
		#if( 1 == length(T0) ) T0 <- rep(T0,nPop)
		#if( 1 == length(TEnd) ) TEnd <- rep(TEnd,nPop)
		if( is.numeric(upperParBounds) ) upperParBounds <- lapply(1:nPop, function(iPop) upperParBounds )
		if( is.numeric(lowerParBounds) ) lowerParBounds <- lapply(1:nPop, function(iPop) lowerParBounds )
		nChains <- dim(x)[3]
		nChainPop <- nChains %/% nPop
		isX <- (0 != length(X))
		isLogDenCompX <- (0 != length(logDenCompX))
		#isTProp <- (0 != length(TProp))
		#if( isTProp && !is.matrix(TProp) ) TProp <- matrix( TProp, nrow=length(TProp), ncol=nPop, dimnames=list(resComp=names(TProp),pop=NULL) )
		pops <- lapply( 1:nPop, function(iPop){
				chainsPopI <- (iPop-1)*nChainPop + 1:nChainPop
				pop <- list( 
					parms=x[,,chainsPopI, drop=FALSE]
					#,T0=T0[iPop]
					#,TEnd=TEnd[iPop]
					,upperParBounds=upperParBounds[[iPop]]
					,lowerParBounds=lowerParBounds[[iPop]]
				)
				if( isX ) pop$X <- X[,chainsPopI , drop=FALSE]
				if( isLogDenCompX ) pop$logDenCompX <- logDenCompX[,chainsPopI ,drop=FALSE]
				#if( isTProp ) pop$TProp <- TProp[,iPop ,drop=TRUE]
				#pop$spaceInd <- spacePop[iPop]
				pop
			})
		#res <- twDEMCBlockInt(pops)
		res <- twDEMCBlockInt(pops=pops,...)
		#.dots <- list(...)
		#do.call( twDEMCBlockInt, c(list(pops=pops),.dots))
		#resc <- concatPops.twDEMCPops(res)	# all have the same length, so allow concatenate population results
	})
attr(twDEMCBlock.array,"ex") <- function(){
	data(twLinreg1)
	attach( twLinreg1 )
	
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	.nPop = 2
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChainPop=4, nPop=.nPop)
	dim(Zinit)
	
	.nGen=100
	#nGen=3
	#mtrace(twDEMC.array)
	#mtrace(.updateIntervalTwDEMCPar)
	#mtrace(twDEMCBlockInt)
	tmp1 <- tmp <-  twDEMCBlock( Zinit, nPop=.nPop
		,dInfos=list(list(fLogDen=logDenGaussian, argsFLogDen=argsFLogDen))
		,nGen=.nGen 
	)
	plot( as.mcmc.list(tmp), smooth=FALSE )
	tmp2 <- tmp <- twDEMCBlock( tmp1, nGen=200 )
	
	str(tmp)
	
	
}

.concatTwDEMCRuns <- function(
	x		##<< initial run (of class twDEMCPops)
	,res	##<< extended run
	, doRecordProposals=FALSE		##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y
){
	##details<< 
	## initial populations are appended by results
	# iPop = length(x$pops)
	for( iPop in 1:length(x$pops) ){
		xPop <- x$pops[[iPop]]
		resPop <- res$pops[[iPop]]
		if( nrow(resPop$parms) < 2){
			# no samples to append
			res$pops[[iPop]] <- xPop	
		}else{
			res$pops[[iPop]]$parms <- abind(xPop$parms, resPop$parms[-1,, ,drop=FALSE], along=1)
			res$pops[[iPop]]$temp <- abind(xPop$temp, resPop$temp[-1, ,drop=FALSE], along=1 )
			res$pops[[iPop]]$pAccept <- abind(xPop$pAccept, resPop$pAccept[-1,, ,drop=FALSE], along=1)
			res$pops[[iPop]]$resLogDen <- abind( xPop$resLogDen, resPop$resLogDen[-1,, ,drop=FALSE], along=1)
			res$pops[[iPop]]$logDen <- abind( xPop$logDen, resPop$logDen[-1,, ,drop=FALSE], along=1)
			nrY <- nrow(resPop$Y)
			if( doRecordProposals || (nrY < 128) ){
				# if new Y has less than 128 rows append previous Y
				res$pops[[iPop]]$Y <- abind( xPop$Y, resPop$Y, along=1)
				nrY <- nrow(res$pops[[iPop]]$Y)	# rows did change
				if( !doRecordProposals && (nrY > 1) )
					res$pops[[iPop]]$Y <- res$pops[[iPop]]$Y[ 1:min(128,nrY),, ,drop=FALSE] 	# cut to 128 cases
			}
		} # if res has samples 
	}
	res
}

setMethodS3("twDEMCBlock","twDEMCPops", function(
		### initialize \code{\link{twDEMCBlockInt}} by former run and append results to former run
		x, ##<< list of class twDEMC, result of \code{\link{twDEMCBlockInt}}
		... ##<< further arguments to \code{\link{twDEMCBlockInt}}
		, TEnd=numeric(0)		##<< numeric vector (nResComp) of end temperatures for each result component, if not given stays at current temperature, if in addition TSpec is specified, it has no effect 
		, doRecordProposals=FALSE		##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y
		, extendRun=TRUE				##<< if set to FALSE then only the new samples are returned
	){
		#twDEMCBlock.twDEMCPops
		##details<< 
		## pops, dInfos, and blocks are reused from x or overwritten by arguments
		nSamplePop <- getNSamples(x)	# used in several cases
		
		.dots <- argsList <- list(...)
		argsList$pops <- x$pops
		argsList$doRecordProposals <- doRecordProposals
		if( 0 == length(.dots$dInfos) )		 
			argsList$dInfos <- x$dInfos
		if( 0 == length(.dots$blocks) )		 
			argsList$blocks <- x$blocks
		#for( iPop in seq_along(x$pops) ){
		#	if( 0 == length(argsList$pops[[iPop]]$T0) )
		#		argsList$pops[[iPop]]$T0 <- 
		#			x$pops[[iPop]]$temp[ nSamplePop[iPop] ]
		#}
		if( 0 == length(argsList$TSpec) ){
			TCurr <- x$pops[[1]]$temp[ nrow(x$pops[[1]]$temp), ]
			if( 0 == length(TEnd) ) TEnd <- TCurr
			if( 1 == length(TEnd) ) TEnd <- rep(TEnd, length(TCurr) )
			argsList$TSpec <- cbind( T0=TCurr, TEnd=TEnd )
		}
		
		if( 0 == length(.dots$controlTwDEMC) )		 
			argsList$controlTwDEMC <- list()
		if( 0 == length(argsList$controlTwDEMC$initialAcceptanceRate) ){
			ar <- abind( lapply( seq_along(x$pops), function(iPop){ 
						adrop( x$pops[[iPop]]$pAccept[ nSamplePop[iPop],, ,drop=FALSE ],1)
					}),along=2)					
			argsList$controlTwDEMC$initialAcceptanceRate <- ar
		}
		res <- res0 <- do.call( twDEMCBlockInt, argsList )
		if( extendRun ) res <- .concatTwDEMCRuns(x,res0,doRecordProposals=doRecordProposals)
		res
	})
attr(twDEMCBlock.twDEMCPops,"ex") <- function(){
	data(twdemcEx1) 		# previous run of twDEMCBlock
	class(twdemcEx1)
	twdemcEx1$thin			# thinning interval
	(nGen0 <- getNGen(twdemcEx1))		# number of generations
	
	# extend by 16 generations
	nGen <- 16
	#mtrace(twDEMCBlock.twDEMCPops)
	res <- twDEMCBlock( twdemcEx1, nGen=nGen )
	
	identical( nGen0+nGen, getNGen(res) )
	plot( as.mcmc.list(res), smooth=FALSE )
}

setMethodS3("twDEMCBlock","twDEMC", function(
		### initialize \code{\link{twDEMCBlockInt}} by former run and append results to former run
		x, ##<< list of class twDEMC, result of \code{\link{twDEMCBlockInt}}
		... ##<< further arguments to \code{\link{twDEMCBlockInt}}
	){
		#twDEMC.twDEMC
		.dots <- list(...)
		argsList <- list(x=x$parms)	
		M0 <- nrow(x$rLogDen)
		#extract X, logDenX, and logDenCompX from parms, but only if not given with \dots
		if( is.null(.dots$X) )		#use is.null, because if provided zero length vector, we want to use it 
			argsList$X <- adrop(x$parms[M0,,,drop=FALSE],2)
		if( is.null(.dots$logDenX) ) 
			argsList$logDenX <- x$rLogDen[M0,,drop=TRUE]
		if( is.null(.dots$logDenCompX) ) 
			argsList$logDenCompX <- adrop(x$logDenComp[M0,,,drop=FALSE],2)
		if( is.null(.dots$nPop) ) 
			argsList$nPop <- getNPops.twDEMC(x)
		res <- do.call( twDEMCBlock.array, c(argsList,.dots))
		res$rLogDen[1:M0,] <- x$rLogDen
		res$logDenComp[,1:M0,] <- x$logDenComp
		res$pAccept[1:M0,] <- x$pAccept
		res$temp[1:M0,] <- x$temp
		res
		### parms appended with the further generations of \code{\link{twDEMCBlockInt}} 
	})
#mtrace(twDEMC.twDEMC)
#mtrace(twDEMCBlockInt)
#twUtest(twDEMC,"test.ZinittwDEMC")




