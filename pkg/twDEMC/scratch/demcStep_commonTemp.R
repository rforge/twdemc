#save code before adapting to differential Temp for different data streams.

twDEMCInt <- function(
	### Differential Evolution Markov Chain
	Zinit
	### initial population: a matrix of number of parameters by number of individuals (d x m0 x nChains) see details
	, nGen = 10
	### number of generations, i.e. steps of the chain
	, ...	
	### further arguments passed to fLogDen
	, controlTwDEMC = list()
	### DEMCzsControl parameters influencing the update and the convergens of the chains (see details)	
	, X=NULL
	### initial active population that will be evolved: a matrix of number of parameters by number of chains (d x N)
	###   if null, it will be initialized to the last row of Zinit
	, logDenX=NULL
	### logDensity of initial values, may save time if already calculated
	, fLogDen
	###	\code{function(theta, ...)} calculates a vector of logDensitys corresponding to different data streams of parameter vector theta 
	### or \code{function(theta, logDenCompX, metropolisStepTemp, ...)} to handle first steps in Multi-step Metropolis decision internally. See details.  
	, argsFLogDen=list()
	### further arguments passed to fLogDen
	, logDenCompX=NULL #numeric(0)  #Zinit[FALSE,1,]
	### numeric matrix: number of internal logDen components x chains, see details
	, fLogDenScale=1
	### scalar multiplied to the result of fLogDen 
	###   allows using functions of negative LogDensity (-1) or Gaussian misfit function (-1/2) instead of logDensity
	, debugSequential=FALSE
	### if TRUE apply is used instead of sfApply, for easier debugging
	, remoteDumpfileBasename=NULL
	### the basename of a dumpfile that is created on error on remote process 
	, nPops = 1
	### the number of populations that do not mix: must be a factor of dimension nChains of Zinit
	, fDiscrProp=NULL
	### function applied to proposal, e.g. to round proposals to to discrete possible values function(theta,...)
	, argsFDiscrProp=list()
	### further arguments to fDiscrProp
	,doRecordProposals=FALSE	##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned
){
	##details<<  
	## This is the central method for applying a Differential Evolution Markov Chain.
	## It is invoked usually by \code{twDEMC} (\code{\link{twDEMC.array}} or \code{\link{twDEMC.twDEMC}})
	
	##seealso<< 
	## Further methods deal with \itemize{
	## \item{ Batch invokation of twDEMC and storing intermediate results: \code{\link{twDEMCBatchInt}}  } 
	## \item{ Generating an initial population for twDEMC: \code{\link{initZtwDEMCNormal}}  }
	## \item{ Transforming the results of twDEMC: \code{\link{subChains.twDEMC}}  }
	## \item{ Transforming the parameter space: \code{\link{transOrigPopt.default}}  }
	## \item{ Invoking fLogDen with proposal in a parallel load balanced way: \code{\link{twCalcLogDenPar}}  }
	## }
	## \code{\link{calcDEMCTemp}}
	## \code{\link{logDenGaussian}}
	## \code{\link{.generateXPropThin}}
	## \code{\link{.doDEMCSteps}}
	## \code{\link{.doDEMCStep}}
	
	##details<<  
	##
	## This method is based on Code of ter Braak  ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain 
	## with snooker updater and fewer chains. Statistics and Computing
	## http://dx.doi.org/10.1007/s11222-008-9104-9 .
	
	##    correspondens parameters and symbols in paper
	##    \tabular{cc}{
	##     d$chains       \tab N \cr
	##     d$parms          \tab d \cr
	##     gamma_par     \tab gamma \cr
	##     gamma_snooker \tab gamma \cr
	##     M0            \tab M0 \cr
	##     mZ          	 \tab M 
	##	  }
	
	##details<< \describe{ \item{Initial population: \code{Z}}{  
	## a matrix of number of parameters by number of individuals (d x m0 x Npop) \describe{
	##	 \item{d}{number of dimensions of parameter vector theta}
	##	 \item{m0}{initials : default 10d/Npop (smaller for fewer parameters) all drawn from prior}
	##   \item{nChains}{number of chains}}
	##   alternatively Zinit may be an twDEMC object, then it will be extended by nGen with parameters of former call
	## }}
	
	##details<< \describe{ \item{Several populations: \code{nPops}}{
	## Chains within population are not independent.
	## In order to assess convergence, one must run several independent populations.
	## This is supported by specifying argument nPops. Then the n chains, i.e. dim(Zinit)[3], are grouped into several populations.
	## }}
	
	##details<< \describe{ \item{Detailed control parameters: \code{controlTwDEMC}}{
	##describe<<  
	ctrl = list( 
		F = 2.38, 		##<< related to multiplicative error (F2=F/sqrt(2*Npar), see eps.mult 
		pSnooker= 0.1,	##<< probability of a snooker update (others parallel updates)
		pGamma1 = 0.1,	##<< probability of jumping to state of another chain (different modes)
		epsMult =0.2,	##<< >0 gives d-dimensional variation around gamma. It adds scaled uncorrelated noise to the proposal. Its advantage over eps.add is that its effect scales with the differences of vectors in the population whereas eps.add does not. if the variance of a dimensions is close to 0, eps.mult gives smaller changes. \cr A uniformly distributed error, i.e. F2*runif(1+-epsMult*prop) multiplied to difference vector from parallel update 
		epsAdd = 0,	   	##<< >0 is needed to ensure that all positions in the space can be reached. For targets without gaps, it can set small or even to 0. \cr sd of normally distributed error added to proposal by parallel or snooker update. 
		thin = 4, 	   	##<< thinning interval	 
		T0=1, Tend=1, 	##<< initial and end Temperature to flatten the density surface, for each population
		pAcceptWindowWidth = 100, ##<< number of generations back over which the acceptance rate is calculated
		probUpDir=0.5 	##<< probability of direction between two states of increasing Density (increase during burin may accelerate convergence)
	)  
	## }}
	ctrl[names(controlTwDEMC)] <- controlTwDEMC
	
	##details<< \describe{ \item{Multistep Metropolis decisions: \code{logDenCompX} }{
	## For performance reasons one may decide to do a first Metropolis decision on parameters
	## so that the possibly costly evaluation of the model functions is not needed.
	## Argument logDenCompX holds the current logDens for those components of return of fLogDen that are handled internally in fLogDen.
	## These will not be regarded in the sum over components Metropolis step in twDEMC.
	## Rather the components of the last accepted step are passed to fLogDen as second argument togheter with the current temperature as third argument.
	## \cr 
	## The rownames of the numeric matrix must correspond to components of result of fLogDen.
	## The number of columns must correspond to the number of chains
	## \cr 
	## Alternatively one can just specify a character vector of return components that are handled internally in the objective function.
	## Then a proper matrix of logDenCompX will be initialized to -Inf, i.e. accepting the next step
	## \cr 
	## A proper initial value for logDenCompX that handels component parms and assures that the step is not rejected
	## because of the initial parameters is:
	## \code{logDenCompX=matrix(-Inf, nrow=1, ncol=.nChains, dimnames=list("parms",NULL))}
	## See the code of \code{\link{logDenGaussian}} for an example multistep version of fLogDen
	## }}
	if( 0 < length(logDenCompX)){
		if( is.character(logDenCompX))
			logDenCompX <- matrix(-Inf*fLogDenScale, nrow=length(logDenCompX), ncol=dim(Zinit)[3], dimnames=list(logDenCompX,NULL))
		if( !is.numeric(logDenCompX) || !is.matrix(logDenCompX) || ncol(logDenCompX)!=dim(Zinit)[3] )
			stop("logDenCompX must be a numeric matrix with one column for each chain and row names correponding to a subst of names of result vector of fLogDen")
	}
	if( !hasArg(controlTwDEMC) ) controlTwDEMC <- list()	#??? why does =list() in declaration ot work
	if( !all(names(controlTwDEMC) %in% names(ctrl)) ){
		warning(paste("unknonw entries",paste(names(controlTwDEMC)[which(!(names(controlTwDEMC) %in% names(ctrl)))],collapse=","),"in controlTwDEMC"))
	}
	d <- as.list(structure(dim(Zinit),names=c("parms","gen","chains")))
	if( d$gen < 4 ) stop(paste("to few initial conditions M0=",d$gen,", m0 = max(4,(8*d)%/%d$chains) is a good choice",sep=""))
	if( (d$chains %% nPops) != 0 ) stop(paste("number of chains dim(Zinit)[3]=",d$chains," must be a multiple of number of populations nPops=",nPops))
	
	#preallocate output of parameters, calculated LogDen, indices of accepted steps, and next proposal
	nThinnedGen = max(nGen %/% ctrl$thin,1);	
	nGen = max(nThinnedGen * ctrl$thin,1)  #do not need to move beyond thinning interval
	Nz <- d$gen+(nThinnedGen)
	Z <- array(NA, dim=c(d=d$parms,nStep=Nz,nChain=d$chains), dimnames=dimnames(Zinit) )
	if( is.null(dimnames(Z))) dimnames(Z) = list( parms=paste("p",1:d$parms,sep="_"), steps=NULL, chains=NULL )
	rLogDen = pAccept = matrix(NA,nrow=Nz, ncol=d$chains)
	temp =  matrix(NA, nrow=Nz, ncol=nPops)
	# record of proposals and rLogDen results
	Y <- NULL	#will be initialized on first result of fLogDen (now we do not know number of result components) 
	
	##details<< \describe{ \item{Initial state: \code{X}}{
	## If initial state X is not specified, the last column (generation) of Z is utilized.
	## If in addition to X, logDenX is specified, the fLogDen will not be avaluated for the initial state of the chains.
	## All the results of fLogDen for the initial state must be finite.
	## }}
	Z[,1:d$gen,] <- Zinit
	mZ=d$gen 					#current number of generations
	if( !hasArg(X) | is.null(X) ){ 
		X <- adrop(Zinit[,dim(Zinit)[2],,drop=FALSE],2)	#last row of Zinit
		logDenX <- NULL
	}
	if( !is.numeric(logDenX)){
		#twCalcLogDenPar expects parameters in columns, need transpose
		.logDenCompXPar <- if( (0 == length(logDenCompX))) NULL else t(logDenCompX) #transpose will fail on NULL
		.resLogDenPar <- twCalcLogDenPar(fLogDen=fLogDen, xProp=t(X), logDenCompX=.logDenCompXPar, argsFLogDen=argsFLogDen, fLogDenScale=fLogDenScale
			,debugSequential=debugSequential, remoteDumpfileBasename=remoteDumpfileBasename, ...)
		logDenX <- .resLogDenPar$logDen
		logDenCompX <- if( is.null(.resLogDenPar$logDenComp)) NULL else t(.resLogDenPar$logDenComp)
	}
	rLogDen[d$gen,] <- logfitness_X <-	logDenX	#initial logDensity
	tmp2 <- !is.finite(logfitness_X)
	if( any(tmp2) ) 
		stop(paste("non-finite logDensity of starting values for chains ",paste(which(tmp2),collapse=","),sep=""))
	nChainsPop <- d$chains %/% nPops 
	
	#calculate temperature steps: exponentially decreasing from T0 to Tend
	if( is.null(ctrl$T0)) ctrl$T0=1
	if( length(ctrl$T0) != nPops) ctrl$T0=rep(ctrl$T0[1],nPops)  #provide Temp for each population
	if( length(ctrl$Tend) != nPops) ctrl$Tend=rep(ctrl$Tend[1],nPops)  #provide Temp for each population
	temp[d$gen,] <- ctrl$T0
	TstepFixed = matrix( sapply( 1:nPops, function(iPop){ pmax(1,calcDEMCTemp( ctrl$T0[iPop], ctrl$Tend[iPop], nGen ))}), nrow=nGen, ncol=nPops)
	
	# initialize further parameters to parallel and snooker steps
	ctrl$Npar12  =(d$parms - 1)/2  # factor for Metropolis ratio DE Snooker update
	ctrl$F2 = ctrl$F/sqrt(2*d$parms)
	ctrl$F1 = 1
	ctrl$pIndStep <- 1.5 #independent state about after 1.5 accepted steps
	ctrl$gInd <- 10*d$parms/nChainsPop	#number of independent generations to sample from
	
	#d$parms <- 16; nChainsPop=8; nPops=2; ctrl=list(thin=8); nGen=200
	# number of requried independent generations to choose from (10*d independent states: TerBraak06 TerBraak08) devided by number of chains
	ar <- rep(0.08, nPops)       #assume initial acceptancer ratio of 0.05 for each population
	nGenBack <- pmin(mZ,ceiling(ctrl$gInd * pmax(1,ctrl$pIndStep/(ar * ctrl$thin))))	#number of genrations to select states from for each population, ctrl$gInd is multiplied by number of rows for one step depending on acceptance rate and thinning but at least one  
	
	##details<< \describe{ \item{Acceptance rate}{
	## The acceptance rate is tracked for each chain across ctrl$pAcceptWindowWidth generations.
	## \cr If acceptance rates drops to low values, this might be because of bad localization
	## ,i.e the requirement to use states from more distant past.
	## In order to improve localization, less parameters or more chains per population are required.
	## }}
	# acceptWindow holds the number of accepted steps for each thinning interval for each chain
	# if its boundaries are filled up, the last nGenBack states are copied to the first interval
	aThinWinW <- ctrl$pAcceptWindowWidth %/% ctrl$thin	# number of thinning intervals over which to calculate average acceptance rate
	ctrl$pAcceptWindowWidth <- aThinWinW * ctrl$thin
	acceptWindow <- matrix( NA, nrow=2*aThinWinW, ncol=d$chains )	#record number of accepted steps in thinning interval
	
	# arguments to .doDEMCStep that do not change within thinning interval
	# tmp <- try( list(...) ); if( inherits(tmp, "try-error")) recover()
	
	.dots <- list(...)
	if( 0<length(.dots))
		if( any(is.null(names(.dots)))) stop("twDEMCInt: encountered \\dots without name. Check ',' after last list entries. Check for using '<-' instead of '='")
	argsDEMCStep = list(
		fDiscrProp=fDiscrProp, argsFDiscrProp=argsFDiscrProp,
		fLogDen=fLogDen,
		fLogDenScale=fLogDenScale,
		argsFLogDen=c(argsFLogDen,.dots)
	)
	argsDEMCStepWrapper <- list( 
		remoteFun=.doDEMCStep,
		argsDEMCStep=argsDEMCStep
	)
	argsDEMCStepWrapper$remoteDumpfileBasename<-remoteDumpfileBasename	#if null delete
	if( !debugSequential & sfParallel() )
		sfExport("argsDEMCStepWrapper")	#evaluated in remote process, only passed one time
	#else( assign("argsDEMCStepWrapper", argsDEMCStepWrapper, pos=1))	#export to current 
	
	# arguments to to demc that change between thinning intervals
	# substract logDen components of logDenCompX from logDenXExt
	chainState <- list(
		X = X
		,logDenCompX = logDenCompX
		,logDenX = logDenX
		,logDenXExt = logDenX -  {if( 0<length(logDenCompX) && all(is.finite(logDenCompX)) ) 
				max(0,sum(logDenCompX)*fLogDenScale) else 0} 	
	)
	
	for( iThin0 in (0:(nThinnedGen-1)) ){
		# calculate proposed steps (differences not destinations) within next thinning interval
		genPropRes <- .generateXPropThin(nPops=nPops, Z=Z,rLogDen=rLogDen,mZ=mZ,ctrl=ctrl,nGenBack=nGenBack)
		
		iGen = iThin0*ctrl$thin+(1:ctrl$thin)
		tempThin = t(TstepFixed[iGen,rep(1:nPops,each=nChainsPop),drop=FALSE])	#chains must be first dimension in order to acces by temp[i]
		#tmp.argsDEMCStep <- if( !debugSequential ) as.name("argsDEMCStep") else argsDEMCStep 
		#resDo <- do.call(.doDEMCSteps, c(chainState[c("X", "logDenCompX", "logDenX")], genPropRes, list(temp=tempThin, argsDEMCStep=tmp.argsDEMCStep, debugSequential=debugSequential)), quote=TRUE )
		#--- debugging
		#debugSequential=TRUE
		#.doDEMCSteps <- twDEMC:::.doDEMCSteps
		#mtrace(.doDEMCSteps)
		# tmpf <- argsDEMCStepWrapper$remoteFun; mtrace(tmpf); argsDEMCStepWrapper$remoteFun<-tmpf
		tmp.remoteFunArgs <- if( !debugSequential & sfParallel() ) as.name("argsDEMCStepWrapper") else argsDEMCStepWrapper 	# eval.parent fails for argsDEMCStepWrapper within sequential function
		resDo <- do.call(.doDEMCSteps, c(chainState[c("X", "logDenCompX", "logDenX", "logDenXExt")], genPropRes, list(temp=tempThin, remoteFunArgs=tmp.remoteFunArgs, debugSequential=debugSequential, doRecordProposals=doRecordProposals)), quote=TRUE )
		chainState <- resDo[ names(chainState) ]
		if( doRecordProposals ){
			if(is.null(Y)){
				.dim <- dim(resDo$Y)
				.dim[2] <- nGen
				Y <- array( double(nGen), dim=.dim, dimnames=c(dimnames(resDo$Y)[1],dimnames(Zinit)[2:3]) )
			} 
			Y[,iGen,] <- resDo$Y
		}
		
		#XXTODO: adapt acceptPos 
		#row in acceptance Window to record acceptance, if exceeds window, copy second part to first (rewind)
		acceptPos0 <- aThinWinW + (iThin0 %% aThinWinW)
		if( acceptPos0 == aThinWinW ){	#interval exeeded or new, copy second half to first half
			acceptWindow[ 1:aThinWinW, ] <- acceptWindow[ aThinWinW+(1:aThinWinW), ]
			acceptWindow[ aThinWinW+(1:aThinWinW), ] <- NA			
		}
		acceptWindow[ acceptPos0+1, ] <- resDo$accept
		
		#after thinning interval add current chain states to past record and adapt variables that depend on past states
		#if (!(iGen%%ctrl$thin) ){   #if modulo is zero
		mZ = mZ+1
		#store result
		Z[,mZ,] <- chainState$X	
		rLogDen[mZ,] <- chainState$logDenX
		#calculate acceptance rate over last min(ctrl$pAcceptWindowWidth, iGen) generations
		curAcceptRows <- (acceptPos0+1)-(min(aThinWinW-1,iThin0):0)
		pAccept[mZ,] <- colSums(acceptWindow[curAcceptRows,,drop=FALSE])/(length(curAcceptRows)*ctrl$thin) 
		temp[mZ,] <- TstepFixed[iThin0*ctrl$thin+ctrl$thin,]
		
		# update the number of generations to choose from, which depends on acceptance rate
		tmp.chains <- rep(1:nPops, each=nChainsPop)
		ar <- pmax(0.05,tapply( pAccept[mZ,], tmp.chains, mean ))	#mean acceptance rate per population, assume minimum 5% to ensure localization (else may go to 0 -> entire history -> no acceptance) 
		nGenBack <- pmin(mZ,ceiling(ctrl$gInd * pmax(1,ctrl$pIndStep/(ar * ctrl$thin))))	#number of genrations to select states from for each population, ctrl$gInd is multiplied by number of rows for one step depending on acceptance rate and thinning but at least one  
		#}		
		
	} # iThin0 in 0:(nThinnedGen-1)
	#res <- list(parms=Z, rLogDen=rLogDen, thin=ctrl$thin, pAccept=pAccept, call=ccl )
	res <- list(parms=Z, rLogDen=rLogDen, thin=ctrl$thin, pAccept=pAccept, temp=temp, logDenCompX = chainState$logDenCompX, Y=Y )
	class(res) <- c( class(res), "twDEMC" )	#monte carlo proposal list
	res
	### list of class \code{twDEMC} (with \code{nStep = M0+nGen%/%thin}) \describe{
	### \item{parms}{ array (d x nStep x nChain) the initial parameters (last row of M0) and the accepted parameter combinations}
	### \item{rLogDen}{ array (nStep x nChain) the logDen of the accepted parameter combinations}
	### \item{pAccept}{ array (nStep x nChain) the acceptance probability over previous ctrl$pAcceptWindowWidth steps }
	### \item{temp}{ vector (nStep) the Temperature used at the step }
	### \item{logDenCompX}{ array (nInternalResult x nChain) the subset of fLogDen results that are handled internally }
	### \item{thin}{ numeric: thinning interval }
	### \item{Y}{ numeric matrix: if(doRecordProposals) record of all proposals together with results components }
	### }
}

.doDEMCSteps <- function(
	### Perform Metropolis steps within next thinning interval.
	X,				##<< numeric matrix: current location of chains rows: parameters, columns: chains
	logDenCompX,	##<< numeric array: result of fLogDen for current location, rows: result components, columns chains  
	logDenX, 		##<< numeric vector of logDen of current position of chains
	logDenXExt,		##<< numeric vector of logDen exluding components handled internally in fLogDen
	xStep, 			##<< array rows: difference vectors in parameter space, cols: chains, zdim: steps
	rExtra,			##<< numeric matrix: row: chain, col: step within thinning interval
	temp,			##<< numeric matrix: row: chain, col: step within thinning interval
	#argsDEMCStep,	##<< see \code{\link{.doDEMCStep}}
	remoteFunArgs,	
	### see \code{\link[twSnowfall]{sfRemoteWrapper}}
	### Must include entries \itemize{
	### \item{remoteFun=.doDEMCStep}
	### \item{remoteDumpfileBasename} 
	### \item{argsDEMCStep}}
	### Can be a name of a previously exported variable.
	debugSequential=FALSE	##<< see \code{\link[twSnowfall]{sfFArgsApplyDep}}
	,doRecordProposals=FALSE	##<< if TRUE then proposals and results of rLogDen are recorded in result$Y.
){
	##seealso<< 
	## \code{\link{twDEMCInt}}
	## \code{\link{.doDEMCStep}}
	
	##details<< 
	## The step must be the last dimension in all arguments in order to make use of dependence step 
	## in load balanced execution.
	xStepStacked <- do.call( rbind, lapply(1:(dim(xStep)[3]),function(iStep){t(adrop(xStep[,,iStep,drop=FALSE],3))}) )
	#all chains of first step, all chains of second step, ...
	d <- as.list(structure(dim(xStep),names=c("parms","chains","steps"))) 
	iGenT <- (1:d$steps)
	nCases = d$chains * d$steps
	if( !(nCases == length(rExtra)) || !(dim(rExtra)==dim(temp)) )
		stop("number of cases in steps must correspond to length of rExtra and length of temp")
	#F_APPLY <- .doDEMCStep	#the function executed on the nodes: one metropolis step
	F_APPLY <- sfRemoteWrapper	#the function executed on the nodes: one metropolis step
	#fArgsSame <- list( remoteFunArgs=as.name("argsDEMCStepWrapper") )	#exported, includes remoteFun=.doDEMCStep, remoteDumpfileBasename and argsDEMCStep
	F_ARGS <- function(i,prevRes){ args<-c(	
			prevRes[c("x", "logDenCompX", "logDenX","logDenXExt")],
			list( step=xStepStacked[i,,drop=TRUE], rExtra=rExtra[i], temp=temp[i]),
			#list( argsDEMCStep=argsDEMCStep )
			list( remoteFunArgs=remoteFunArgs )
		)}	
	#.res0 <- lapply(1:nrow(X),function(row){X[row,]})
	res0 <- lapply(1:d$chains,function(iChain){list(
				x=X[,iChain,drop=TRUE],
				logDenCompX=if(0<length(logDenCompX)) logDenCompX[,iChain,drop=TRUE] else logDenCompX
				,logDenX=logDenX[iChain]
				,logDenXExt=logDenXExt[iChain]
			)})
	res <- sfFArgsApplyDep( nCases, F_ARGS, F_APPLY, res0, debugSequential=debugSequential)
	#all chains of first step, all chains of second step, ...
	#iChains <- (d$steps-1)*d$chains+(1:d$chains)	#index of the chains of the last step
	#modify in place, so that dimnames etc are preserved
	endChain0 <- (d$steps-1)*d$chains	#index before last step of all chains
	#s0 <- list(X=X, logDenCompX=logDenCompX, logDenX=logDenX)#save former state
	boResFLogDenX <- { .nc <- ncol(logDenCompX); ( !is.null(.nc) && (.nc > 0) ) }	#if number of columns > 0
	#logDenCompComp
	for( iChain in 1:d$chains ){
		resChain = res[[endChain0+iChain]]
		X[,iChain] <- resChain$x
		if(	boResFLogDenX ) logDenCompX[,iChain] <- resChain$logDenCompX
		logDenX[iChain] <- resChain$logDenX
		logDenXExt[iChain] <- resChain$logDenXExt
	}
	Y <- NULL
	if( doRecordProposals){ 
		parNames <- rownames(X)
		resCompNames <- names(res[[1]]$logDenCompComp)
		tmp <- c("rLogDen",parNames,"accepted",resCompNames) 
		Y <- array( double(length(tmp)*d$steps*d$chains), dim=c(d=length(tmp),nStep=d$steps,nChain=d$chains),	dimnames=list(comp=tmp,steps=NULL,chains=NULL) )			
		for(iStep in 1:d$steps){
			chain0 <- (iStep-1)*d$chains
			for( iChain in 1:d$chains ){
				i <- chain0+iChain
				Y["accepted",iStep,iChain] <- res[[i]]$accepted
				Y[parNames,iStep,iChain] <- res[[i]]$xProp
				Y[resCompNames,iStep,iChain] <- resComp <- res[[i]]$logDenCompComp
				Y["rLogDen",iStep,iChain] <- sum(resComp)
			}
		}
	}
	acceptedM <- matrix( sapply(res,function(resChain){ resChain$accepted}), byrow=TRUE, ncol=d$chains )	#rows: steps, cols: chains
	accepted <- colSums(acceptedM)
	
	resDo <- list(	X=X, logDenCompX=logDenCompX, logDenX=logDenX, logDenXExt=logDenXExt, accepted=accepted, Y=Y )
	### list with components \describe{
	### \item{X}{matrix current position, column for each chain}
	### \item{logDenCompX}{matrix: internal result components of fLogDen current position, column for each chain}
	### \item{logDenX}{vector current logDen of chains}
	### \item{logDenX}{vector current logDen of chains, excluding components of logDen that are handled within fLogDen}
	### \item{accepted}{numerical vector, number of accepted steps for each chain}
	### \item{Y}{numerical matrix (steps x 2+params+result): accepted, rLogDen, parms, and all fLogDen result components for each proposal }
}


.doDEMCStep <- function( 
	### Perfrom one DEMC step, function to be  alled in remote process.
	x				##<< numeric vector: current state
	, logDenCompX	##<< numeric vector: currnet states compoenents of fLogDen results
	, logDenX		##<< numeric scalar
	, logDenXExt	##<< numeric scalar
	,step			##<< numeric vector of step in X space
	,rExtra 		##<< corrector for rLogDen due to selection of step 
	,temp			##<< temperature of the current step and population
	,argsDEMCStep	
### list with components \describe{
### \item{fDiscrProp,argsFDiscrProp}{function and additional arguments applied to xProp, e.g. to round it to discrete values}
### \item{argsFLogDen, fLogDenScale}{additional arguments to fLogDen and scalar factor applied to result of fLogDen}
### \item{}{}
### }
){
	#.doDEMCStep
	##seealso<< 
	## \code{\link{twDEMCInt}}
	## \code{\link{.doDEMCSteps}}
	boResFLogDenX <- (length(logDenCompX) > 0)
	with( argsDEMCStep, {
			#attach( argsDEMCStep )
			accepted<-FALSE
			xProp = x + step
			if( is.function(fDiscrProp)) xProp = do.call(fDiscrProp,xProp,argsFDiscrProp, quote=TRUE)
			res <- if(boResFLogDenX)
					do.call( fLogDen, c(list(xProp, logDenCompX, temp), argsFLogDen) )	# evaluate logDen
				else
					do.call( fLogDen, c(list(xProp), argsFLogDen) )	# evaluate logDen
			##details<< 
			## if logDenCompX is given, then these components of result of fLogDen are handled
			## internally. Hence, for Metropolis step, constrain res to remaining components
			resOut <- if(boResFLogDenX){
					posLogDenInt <- match(names(logDenCompX), names(res) )	#side effect: used below
					if( any(is.na(posLogDenInt)))
						stop(paste("not all names of logDenCompX in return of fLogDen: ",paste(names(res),collapse=",")))
					res[-posLogDenInt]
				}else
					res
			logDenProp <- sum(res)*fLogDenScale			# sum over all components
			logDenPropExt <- sum(resOut)*fLogDenScale	# sum only over those components that are not handeled with of-internal metroplis decisions
			if( is.finite(logDenPropExt) ){
				#Metropolis step
				logr = (logDenPropExt+rExtra - logDenXExt) / temp
				if ( is.finite(logr) & (logr) > log(runif(1)) ){
					accepted <- TRUE
					x <- xProp
					if(boResFLogDenX) logDenCompX <- res[posLogDenInt]
					logDenX <- logDenProp
					logDenXExt <- logDenPropExt
				}
			}
			list(accepted=accepted,x=x,logDenCompX=logDenCompX,logDenX=logDenX, logDenXExt=logDenXExt,logDenCompComp=res,xProp=xProp)
		})
	#detach( argsDEMCStep ); list(accepted=accepted,x=x,logDenCompX=logDenCompX,logDenX=logDenX)
	### list with components \describe{
	### \item{accepted}{boolean scalar: if step was accepted}
	### \item{x}{numeric vector: current position in parameter space}
	### \item{logDenCompX}{numeric vector: result components of fLogDen for current position that are handeled inside fLogDen}
	### \item{logDenX}{numeric scalar: logDen of current position}}
	### \item{logDenXExt}{numeric scalar: logDen of current position, exlcuding the components that are handled inside fLogDen.}}
	### \item{xProp}{numeric vector: proposal}
	### \item{logDenCompComp}{numeric vector: result components of fLogDen for current position}
}