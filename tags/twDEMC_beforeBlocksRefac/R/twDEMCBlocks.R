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
			, compPos=1:nrow(Zinit)	##<< index or names of the parameter components to be updated
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
			, resCompNames			##<< character vector of result component names
			, intResCompNames=character(0) 
				### character vector: names of results components of fLogDen 
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
	if( any( !is.function( sapply(blocks, "[[", "fLogDen")) ) ) stop("all blocks need function fLogDen")
	if( any( !is.function( sapply(blocks, "[[", "fUpdate")) ) ) stop("all blocks need function fLogDen")
	d <- as.list(structure(dim(Zinit),names=c("parms","gen","chains")))
	if( sort(sapply(blocks,"compPos")) != 1:d$parms ) stop("union of all blocks must equal rows of Zinit")
	
	ctrl = list( 
		thin = 4 	   	##<< thinning interval	 
		,TFix = vector("numeric",0)		##<< named numeric vector: result components for which temperature shoudl not change		
		,F = 2.38 		##<< related to multiplicative error (F2=F/sqrt(2*Npar), see eps.mult 
		,pSnooker= 0.1,	##<< probability of a snooker update (others parallel updates)
		pGamma1 = 0.1,	##<< probability of jumping to state of another chain (different modes)
		epsMult =0.2,	##<< >0 gives d-dimensional variation around gamma. It adds scaled uncorrelated noise to the proposal. Its advantage over eps.add is that its effect scales with the differences of vectors in the population whereas eps.add does not. if the variance of a dimensions is close to 0, eps.mult gives smaller changes. \cr A uniformly distributed error, i.e. F2*runif(1+-epsMult*prop) multiplied to difference vector from parallel update 
		epsAdd = 0,	   	##<< >0 is needed to ensure that all positions in the space can be reached. For targets without gaps, it can set small or even to 0. \cr sd of normally distributed error added to proposal by parallel or snooker update. 
		pAcceptWindowWidth = 256, ##<< number of generations back over which the acceptance rate is calculated
		probUpDir=0.5 	##<< probability of direction between two states of increasing Density (increase during burin may accelerate convergence)
		,initialAcceptanceRate=rep(0.25,nChains)	##<< numeric vector (nChains) initially assumed acceptance rate. Used to calculate the number of generations backwards to sample from
		,DRgamma=0		##<< factor for reducing step length [0..1) in delayed rejection step, 0 means no DR step
		,minPCompAcceptTempDecr=0.15  ##<< if acceptance rate drops below minPCompAcceptTempDecr+0.02 times this level, employ delayed rejection (DR)
		,pIndStep = 1.5 ##<< independent state about after on average about those number of 1.5 accepted steps
		,nPastGen = 10  ##<< factor for determining the number of recent past states to sample during burnin. It is multiplied by the number of parameters. Past generations are calculated by deviding by the number of chains per population 
	)  ##end<< 
	## }}
	ctrl[names(controlTwDEMC)] <- controlTwDEMC
	
	# dimensions of Zinit
	ZinitPops <- lapply(pop,"[[","Zinit")
	parNames <- rownames(ZinitPops[[1]])
	nParm <- nrow(ZinitPops[[1]])
	nPop <- length(pops)
	nChainsPop <- dim(ZinitPops[[1]])[3]
	nChains <- nChainsPop*nPop 
	M0Pops <- sapply( ZinitPops, ncol )
	iPops <- 1:nPop
	
	# dimensions of blocks
	nBlock <- length(blocks)
	iBlocks <- 1:length(blocks)
	compPosBlock <- lapply( blocks, function(block){ 
			compPosBi <- block$compPos
			if( is.character(compPosBi)) match(compPosBi,parNames ) else compPosBi 
		})
	#XX think about invoking block$fLogDen to geth resCompNames instead of argument 
	resCompNamesOrig <- lapply( blocks, "[[", "resCompNames" )
	resCompNamesUnique <- lapply( iBlocks, function(iBlock){ 
			.getRescompNameBlock(resCompNamesOrig[[iBlock]], iBlock)
		})
	resCompNamesFlat <- do.call(c, resCompNamesUnique)	# concatenate names
	nResComp <- length(resCompNamesFlat)
	resCompNamesPos <- { # list: for each block: position of resCompNames in resCompNamesFlat
		tmp.length <- sapply(resCompNamesUnique, length)
		tmp.start <- cumsum(tmp.length)-1
		tmp.i <- lapply(iBlocks,function(iBlock){ tmp.start[iBlock]+1:tmp.length(iBlock) })
	}
	
	#preallocate output of parameters, calculated LogDen, indices of accepted steps, and next proposal
	nGenPops <- sapply( pops, "[[", "nGen")
	nThinnedGenPops = sapply( nGenPops %/% ctrl$thin, max, 1 )	#number of thinning intervals 
	nGenPops = nThinnedGenPops * ctrl$thin  #number of generations in this batch (do not need to move behind last thinning interval)
	Nz <- M0Popsd+(nThinnedGenPops)			   #number of rows in Z after batch
	ZPops <- lapply( iPops, function(iPop){
			dNames <- if( is.null(dimnames(ZinitPops[[iPop]])) ) dimnames(ZinitPops[[iPop]]) else
					list( parms=paste("p",1:nParm,sep="_"), steps=NULL, chains=NULL )
			array(NA_real_, dim=c(nParm,Nz[iPop],nChainsPop), dimnames=dNames)    
		})
	# for each step (rows)
	# one temperature per populations 
	temp =  matrix(NA_real_, nrow=Nz, ncol=nPops)
	# one rLogDen per chain and per block
	rLogDen = pAccept = array(NA_real_, dim=c(Nz, nBlock, nChains)
		, dimnames=list(steps=NULL, block=NULL, chains=NULL) )
	# several resLogDen per chain (blocks are concatenated)
	resLogDen = array( NA_real_, dim=c(Nz, nResComp, nChains)
		,dimnames=list(steps=NULL, resComp=resCompNamesFlat, chains=NULL)) 
	
	
	
	
	NULL
}
attr(twDEMCBlockInt,"ex") <- function(){
	data(twLinreg1); attach( twLinreg1 ) 
	
	# run with saving a restart file each 50 generations
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,				### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		xval=xval
	)
	argsFLogDenA <- c( argsFLogDen, list(
		thetaPrior= thetaTrue["a"],	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta[1,1,drop=FALSE]	### the inverse of the Covariance of the prior parameter estimates
	))
	argsFLogDenB <- c( argsFLogDen, list(
			thetaPrior= thetaTrue["b"],	### the prior estimate of the parameters
			invCovarTheta = invCovarTheta[2,2,drop=FALSE]	### the inverse of the Covariance of the prior parameter estimates
		))
	
	.nPops=2
	Zinit <- ZinitPops <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=4*.nPops, nPops=.nPops)
	pops <- list(
		pop1 <- list(
			Zinit = ZinitPops[,1:3,1:4,drop=FALSE]	# the first population with less initial conditions
			,nGen=10
			),
		pop2 <- list(
			Zinit = ZinitPops[,1:4,5:8,drop=FALSE]	# the first population with less initial conditions
			,nGen=15
		)
		)
	tmp <- .checkPop(pops[[1]])
	
	res <-  twDEMCBlockInt( Zinit=Zinit, nGen=60,	nPops=.nPops
		, fLogDen=logDenGaussian, argsFLogDen=argsFLogDen
		, nBatch=50	# save each 50 generations
	)
	
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
	if( is.null(pop$Zinit) || !is.numeric(pop$Zinit) || length(dim(pop$Zinit)) != 3 )
		stop(".checkPop: entry Zinit must be a numeric 3D array")
	if( is.null(pop$nGen) ) pop$nGen <- nGenDefault
	if( is.null(pop$T0) ) pop$T0 <- 1
	if( is.null(pop$Tend) ) pop$Tend <- 1
	if( is.null(pop$X) ) pop$X <- pop$Zinit[,ncol(pop$Zinit),]
	### pop with entries nGen, T0, Tend and X defined 
	pop	
}

updateBlockDEMC <- function(
	### updating parameter block by differential evolution 
){
	
}

