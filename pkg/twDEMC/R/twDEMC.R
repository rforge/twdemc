twDEMCInt <- function(
		### Differential Evolution Markov Chain
	Zinit
		### initial population: a matrix of number of parameters by number of individuals (d x m0 x nChains) see details and \code{\link{initZtwDEMCNormal}}.
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
		###	\code{function(theta, ...)} calculates a vector of logDensities 
		### corresponding to different data streams of parameter vector theta 
		### \cr or \code{function(theta, resFLogDenX, metropolisStepTemp, ...)} 
		### to handle first steps in Multi-step Metropolis decision internally. 
		### See details.  
	, argsFLogDen=list()
		### further arguments passed to fLogDen
	, resFLogDenX=NULL #numeric(0)  #Zinit[FALSE,1,]
		### numeric matrix: logDen components x chains, see details
	, intResCompNames=character(0)	
		### character vector: names of results components of fLogDen that are used for internal Metropolis decisions 
	, fLogDenScale=1
		### scalar multiplied to the result of fLogDen 
		### allows using functions of negative LogDensity (-1) or Gaussian misfit function (-1/2) instead of logDensity
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
	, doRecordProposals=FALSE	##<< if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y
){
	##alias<< twDEMC
	
	##details<<  
	## This is the central method for applying a Differential Evolution Markov Chain given a function of 
	## an unnormalized probability density.
	## It is invoked usually by \code{twDEMC} (\code{\link{twDEMC.array}} or \code{\link{twDEMCBatchInt}})
	
	##seealso<< 
	## Further functionality of the twDEMC package deals with \itemize{
	## \item{ Batch invokation of twDEMC and storing intermediate results: \code{\link{twDEMCBatchInt}}  } 
	## \item{ Generating an initial population for twDEMC: \code{\link{initZtwDEMCNormal}}  }
	## \item{ Transforming the results of twDEMC: \code{\link{subChains.twDEMC}}  }
	## \item{ Transforming the parameter space: \code{\link{transOrigPopt.default}}  }
	## \item{ Invoking fLogDen with proposal in a parallel load balanced way: \code{\link{twCalcLogDenPar}}  }
	## \item{ Plotting routines: \code{\link{plotMarginal2D}}  }
	## \item{ Calculating marginal aggregates: \code{\link{marginals1d}}  }
	## \item{ Diagnostics for MCMC samples: \code{\link{checkConvergenceGelman}}  }
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
	
	if( length(nGen) != 1 ) stop("twDEMCInt: nGen must be a scalar.")

	d <- as.list(structure(dim(Zinit),names=c("parms","gen","chains")))
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
		Tprop=1,		
			### numeric matrix (result components x populations) proportions of Temperature for different data streams
			### Alternatively twDEMCInt can handle vectors (components) which are replicated for populations.
		useMultiT=TRUE,	##<< if TRUE Temperature for different data streams are scaled by logDen of accepted state
		TFix=numeric(0),	##<< named numeric vector of Temperatures of fLogDen components, whose Temperature is fixed
		pAcceptWindowWidth = 256, ##<< number of generations back over which the acceptance rate is calculated
		probUpDir=0.5 	##<< probability of direction between two states of increasing Density (increase during burin may accelerate convergence)
		,initialAcceptanceRate=rep(0.25,nPops)	##<< numeric vector (nPops) initially assumed acceptance rate. Used to calculate the number of generations backwards to sample from
		,DRgamma=0		##<< factor for reducing step length [0..1) in delayed rejection step, 0 means no DR step
		,minPCompAcceptTempDecr=0.2  ##<< if acceptance rate drops below minPCompAcceptTempDecr+0.02 times this level, employ delayed rejection (DR)
	)  ##end<< 
	## }}
	ctrl[names(controlTwDEMC)] <- controlTwDEMC
	
	##details<< \describe{ \item{Multistep Metropolis decisions: \code{intResCompNames} }{
	## For performance reasons one may decide to do a Metropolis decision on parameters
	## before evaluating the model in order to save the possibly costly model evaluation.
	## \cr
	## Argument \code{intResCompNames} specifies a character vector of return components that are handled internally in the objective function.
	## The first currently accepted value provided to next 
	## evaluation of the objective function is taken from \code{resFLogDenX}. If this 
	## argument is NULL, then a proper matrix of \code{resFLogDenX} will be initialized to -Inf, i.e. accepting the next step
	## \cr 
	## See the code of \code{\link{logDenGaussian}} for an example multistep version of fLogDen
	## }}
	.dots <- list(...)
	if( any(""==names(.dots)) || length(names(.dots))!=length(.dots) )
		stop("twCalcLogDenInt: encountered unnamed argument in ... Check for <- and ,, in list()")
	
	tmp.f <- function(){	#XX TODO: for backward compatibitiy  set rseFLogDenX and internalNames
		if( 0 < length(resFLogDenX)){
			if( is.character(resFLogDenX))
				resFLogDenX <- matrix(-Inf*fLogDenScale, nrow=length(resFLogDenX), ncol=dim(Zinit)[3], dimnames=list(resFLogDenX,NULL))
			if( !is.numeric(resFLogDenX) || !is.matrix(resFLogDenX) || ncol(resFLogDenX)!=dim(Zinit)[3] )
				stop("resFLogDenX must be a numeric matrix with one column for each chain and row names correponding to a subst of names of result vector of fLogDen")
		}
	}
	if( !hasArg(controlTwDEMC) ) controlTwDEMC <- list()	#??? why does =list() in declaration ot work
	if( !all(names(controlTwDEMC) %in% names(ctrl)) ){
		warning(paste("unknonw entries",paste(names(controlTwDEMC)[which(!(names(controlTwDEMC) %in% names(ctrl)))],collapse=","),"in controlTwDEMC"))
	}
	if( d$gen < 4 ) stop(paste("to few initial conditions M0=",d$gen,", m0 = max(4,(8*d)%/%d$chains) is a good choice",sep=""))
	if( (d$chains %% nPops) != 0 ) stop(paste("number of chains dim(Zinit)[3]=",d$chains," must be a multiple of number of populations nPops=",nPops))

	#preallocate output of parameters, calculated LogDen, indices of accepted steps, and next proposal
	nThinnedGen = max(nGen %/% ctrl$thin,1);	#number of thinning intervals in this batch
	nGen = max(nThinnedGen * ctrl$thin,1)  #number of generations in this batch (do not need to move behind last thinning interval)
	Nz <- d$gen+(nThinnedGen)			   #number of rows in Z after batch
	Z <- array(NA, dim=c(d=d$parms,nStep=Nz,nChain=d$chains), dimnames=dimnames(Zinit) )
	if( is.null(dimnames(Z))) dimnames(Z) = list( parms=paste("p",1:d$parms,sep="_"), steps=NULL, chains=NULL )
	rLogDen = pAccept = matrix(NA,nrow=Nz, ncol=d$chains)
	temp =  matrix(NA, nrow=Nz, ncol=nPops)
	
	##details<< \describe{ \item{Initial state: \code{X}}{
	## If initial state X is not specified, the last column (generation) of Z is utilized.
	## If in addition to X, logDenX is specified, the fLogDen will not be avaluated for the initial state of the chains.
	## All the results of fLogDen for the initial state must be finite.
	## }}
	Z[,1:d$gen,] <- Zinit
	mZ=d$gen 					#current number of generations
	if( !hasArg(X) | 0 == length(X) ){ 
		# assume current state X of each chain is last row of Zinit, then reset logDen
		X <- adrop(Zinit[,dim(Zinit)[2],,drop=FALSE],2)	#last row of Zinit
		logDenX <- numeric(0)
		resFLogDenX <- numeric(0)
	}
	if( 0==length(rownames(X)) )
		warning("twDEMCInt: no rownames of initial value matrix X.")
	if( (0==length(logDenX)) || (0==length(resFLogDenX)) || !is.finite(logDenX) || any(!is.finite(resFLogDenX))){
		# calculate unnormalized density for current state if not provided yet 
		#twCalcLogDenPar expects parameters in columns, need transpose
		.resFLogDenXPar <- if( (0 == length(resFLogDenX))) numeric(0) else t(resFLogDenX) #transpose will fail on NULL
		.resLogDenPar <- twCalcLogDenPar(fLogDen=fLogDen, xProp=t(X), resFLogDenX=.resFLogDenXPar, intResCompNames=intResCompNames, argsFLogDen=argsFLogDen, fLogDenScale=fLogDenScale
			,debugSequential=debugSequential, remoteDumpfileBasename=remoteDumpfileBasename, ...)
		logDenX <- .resLogDenPar$logDen
		resFLogDenX <- t(.resLogDenPar$resFLogDen)
		rownames(resFLogDenX) <- .getResFLogDenNames(resFLogDenX)
	}else{
		# invoke fLogDen once to check for consistency of result components
		if( 0 < length(intResCompNames)){
			.logDenAcc <- resFLogDenX[intResCompNames,1]
			.temp <- structure(rep(1.0,length(.logDenAcc)),names=names(.logDenAcc))
			tmp <- do.call( fLogDen, c( list(X[,1], .logDenAcc, .temp), argsFLogDen, list(...)) )
		}else
			tmp <- do.call( fLogDen, c( list(X[,1]),argsFLogDen, list(...)) )
		if( any(!is.finite(tmp)))
			stop(paste("twDEMCInt: non-finite logDensity of starting value for chain 1",sep=""))
		if( !identical(.getResFLogDenNames(tmp), .getResFLogDenNames(resFLogDenX)))
			stop(paste("twDEMCInt: non-matching names of resFLogDenX and result of fLogDen",sep=""))
	}
	rLogDen[d$gen,] <- logfitness_X <-	logDenX	#initial logDensity
	tmp2 <- !is.finite(logfitness_X)
	if( any(tmp2) ) 
		stop(paste("twDEMCInt: non-finite logDensity of starting values for chains ",paste(which(tmp2),collapse=","),sep=""))
	nChainsPop <- d$chains %/% nPops 
	
	# record of proposals and fLogDen results, rows c(accepted, parNames, resCompNames, rLogDen)
	# if doRecordProposals is FALSE record only the thinning intervals for the last 128 generations
	# +1 for the initial state
	nThinLast <- if(doRecordProposals) nThinnedGen else min(nThinnedGen, ceiling(128/ctrl$thin))
	nThinOmitRecord = nThinnedGen-nThinLast	#the Thinning intervals with no recording of outputs
	nGenOmitRecord = nThinOmitRecord*ctrl$thin
	nGenY <- nThinLast*ctrl$thin+1
	Y <- array( double((d$parms+nrow(resFLogDenX)+2)*nGenY*d$chains), dim=c( d$parms+nrow(resFLogDenX)+2, nGenY, d$chains)
			 , dimnames=c(list(vars=c("rLogDen",rownames(X),"accepted",rownames(resFLogDenX))),dimnames(Zinit)[2:3]) )
	
	#calculate temperature steps: exponentially decreasing from T0 to Tend
	if( is.null(ctrl$T0)) ctrl$T0=1
	if( length(ctrl$T0) == 1) ctrl$T0=rep(ctrl$T0[1],nPops)  #provide Temp for each population
	if( length(ctrl$Tend) == 1) ctrl$Tend=rep(ctrl$Tend[1],nPops)  #provide Temp for each population
#cat("twDEMCInt: ctrl$Tend=",paste(round(ctrl$Tend),collapse=","),"\n",sep="")	
	temp[d$gen,] <- ctrl$T0
	TstepFixed = matrix( sapply( 1:nPops, function(iPop){ pmax(1,calcDEMCTemp( ctrl$T0[iPop], ctrl$Tend[iPop], nGen ))}), nrow=nGen, ncol=nPops)
#cat("twDEMCInt: TstepFixed[1,]=",paste(round(TstepFixed[1,]),collapse=","),"\n",sep="")	
#cat("twDEMCInt: TstepFixed[nGen,]=",paste(round(TstepFixed[nGen,]),collapse=","),"\n",sep="")	

	# initialize further parameters to parallel and snooker steps
	ctrl$Npar12  =(d$parms - 1)/2  # factor for Metropolis ratio DE Snooker update
	ctrl$F2 = ctrl$F/sqrt(2*d$parms)
	ctrl$F1 = 1
	ctrl$pIndStep <- 1.5 #independent state about after 1.5 accepted steps
	ctrl$gInd <- 10*d$parms/nChainsPop	#number of independent generations to sample from
	
	#d$parms <- 16; nChainsPop=8; nPops=2; ctrl=list(thin=8); nGen=200
	# number of requried independent generations to choose from (10*d independent states: TerBraak06 TerBraak08) devided by number of chains
	ar <- ctrl$initialAcceptanceRate 
	if( length(ar) != nPops)
		stop("ctrl$initialAcceptanceRate must be of length nPops")
	nGenBack <- pmin(mZ,ceiling(ctrl$gInd * pmax(1,ctrl$pIndStep/(ar * ctrl$thin))))	#number of genrations to select states from for each population, ctrl$gInd is multiplied by number of rows for one step depending on acceptance rate and thinning but at least one  
	
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
	acceptWindow <- matrix( rep(ctrl$thin*ar, each=2*aThinWinW*nChainsPop), nrow=2*aThinWinW, ncol=d$chains )	#record number of accepted steps in thinning interval
	
	# arguments to .doDEMCStep that do not change within thinning interval
	# tmp <- try( list(...) ); if( inherits(tmp, "try-error")) recover()

	#match the positions in resFLogDenX that are handled internally
	posLogDenInt <- match(intResCompNames, .getResFLogDenNames(resFLogDenX) )	
	if( any(is.na(posLogDenInt)))
		stop(paste("not all names of intResCompNames (",paste(intResCompNames,collapse=","),") in return of fLogDen: ",paste(rownames(resFLogDenX),collapse=",")))
	
	#proportions of temperature numeric matrix (comp x populations) for temperature varying across data streams
	Tprop=matrix(1, nrow=nrow(resFLogDenX), ncol=nPops, dimnames=dimnames(resFLogDenX))
	if( length(ctrl$Tprop) > 1){
		if( !is.matrix(ctrl$Tprop) )
			ctrl$Tprop <- matrix( ctrl$Tprop, nrow=nrow(resFLogDenX), ncol=nPops, dimnames=dimnames(resFLogDenX) )
		Tprop <- ctrl$Tprop[ .getResFLogDenNames(resFLogDenX), ,drop=FALSE]
		if( nrow(Tprop) != nrow(resFLogDenX))
			stop("ctrl$Tprop must have a named entry for each component of fLogDen")
		Tprop[] = apply(Tprop,2,function(Tprop){Tprop / max(Tprop)}) #scale to maximum 1 per population
	}

	argsDEMCStep = list(
		fDiscrProp=fDiscrProp, argsFDiscrProp=argsFDiscrProp,
		fLogDen=fLogDen,
		fLogDenScale=fLogDenScale,
		argsFLogDen=c(argsFLogDen,.dots)
		,posLogDenInt=posLogDenInt	
		,useMultiT=ctrl$useMultiT
		,DRgamma=ctrl$DRgamma
		,minPCompAcceptTempDecr=ctrl$minPCompAcceptTempDecr
		,Tprop=Tprop
		,TFix=ctrl$TFix
		,posTFix=match(names(ctrl$TFix),.getResFLogDenNames(resFLogDenX))
	)
	argsDEMCStepWrapper <- list( 
		remoteFun=.doDEMCStep,
		argsDEMCStep=argsDEMCStep	# doDEMCStep 
	)
	argsDEMCStepWrapper$remoteDumpfileBasename<-remoteDumpfileBasename	#if null delete
	if( !debugSequential & sfParallel() )
		sfExport("argsDEMCStepWrapper")	#evaluated in remote process, only passed one time
	#else( assign("argsDEMCStepWrapper", argsDEMCStepWrapper, pos=1))	#export to current 
	tmp.remoteFunArgs <- if( !debugSequential & sfParallel() ) as.name("argsDEMCStepWrapper") else argsDEMCStepWrapper 	# eval.parent fails for argsDEMCStepWrapper within sequential function

	# arguments to to demc that change between thinning intervals
	# substract logDen components of resFLogDenX from logDenXExt
	chainState <- list(
		X = X
		,resFLogDenX = resFLogDenX
		,logDenX = logDenX
	)
	
	for( iThin0 in (0:(nThinnedGen-1)) ){
		# calculate proposed steps (differences not destinations) within next thinning interval
		genPropRes <- .generateXPropThin(nPops=nPops, Z=Z,rLogDen=rLogDen,mZ=mZ,ctrl=ctrl,nGenBack=nGenBack)

		iGen = iThin0*ctrl$thin+(1:ctrl$thin)
		tempThin = t(TstepFixed[iGen,rep(1:nPops,each=nChainsPop),drop=FALSE])	#chains must be first dimension in order to acces by temp[i]
		if( iThin0 == nThinOmitRecord ){
			# record the initial state of the proposals record
			Y["rLogDen",1,] <- chainState$logDenX
			Y[rownames(X),1,] <- chainState$X
			Y["accepted",1,] <- 1	#TRUE
			Y[rownames(resFLogDenX),1,] <- chainState$resFLogDenX
		}
		boRecordProposalsIThin = (iThin0 >= nThinOmitRecord)	# only record the last proposals after nThinOmitRecord
		#tmp.argsDEMCStep <- if( !debugSequential ) as.name("argsDEMCStep") else argsDEMCStep 
		#resDo <- do.call(.doDEMCSteps, c(chainState[c("X", "resFLogDenX", "logDenX")], genPropRes, list(temp=tempThin, argsDEMCStep=tmp.argsDEMCStep, debugSequential=debugSequential)), quote=TRUE )
		#--- debugging
		#debugSequential=TRUE
		#.doDEMCSteps <- twDEMC:::.doDEMCSteps
		#mtrace(.doDEMCSteps)
		# tmpf <- argsDEMCStepWrapper$remoteFun; mtrace(tmpf); argsDEMCStepWrapper$remoteFun<-tmpf
		# tmp.remoteFunArgs <- if( !debugSequential & sfParallel() ) as.name("argsDEMCStepWrapper") else argsDEMCStepWrapper 	# eval.parent fails for argsDEMCStepWrapper within sequential function
		resDo <- do.call(.doDEMCSteps, c(chainState[c("X", "resFLogDenX", "logDenX")], genPropRes, list(temp=tempThin
				, nPops=nPops
				, pAccept=ar
				, remoteFunArgs=tmp.remoteFunArgs
				, debugSequential=debugSequential
				, doRecordProposals=boRecordProposalsIThin
			))
			, quote=TRUE )
		chainState <- resDo[ names(chainState) ]
		
		#row in acceptance Window to record acceptance, if exceeds window, copy second part to first (rewind)
		acceptPos0 <- aThinWinW + (iThin0 %% aThinWinW)
		if( acceptPos0 == aThinWinW ){	#interval exeeded or new, copy second half to first half
			acceptWindow[ 1:aThinWinW, ] <- acceptWindow[ aThinWinW+(1:aThinWinW), ]
			acceptWindow[ aThinWinW+(1:aThinWinW), ] <- NA			
		}
		acceptWindow[ acceptPos0+1, ] <- resDo$accept

		#record results of thinning interval
		mZ = mZ+1
		#store result
		Z[,mZ,] <- chainState$X	
		rLogDen[mZ,] <- chainState$logDenX
		#calculate acceptance rate over last min(ctrl$pAcceptWindowWidth, max(iGen,4) ) generations
		#the first 5 thinning generations ctrl$intialAcceptanceRate in part determines the generations to go back
		#curAcceptRows <- (acceptPos0+1)-((aThinWinW-1):0)	# initial phase of batch strongly determined by prescribed assumed initial acceptance rate
		#curAcceptRows <- (acceptPos0+1)-(min(aThinWinW-1,max(iThin0,4)):0) # initial phase strongly determined by random first proposals
		curAcceptRows <- (acceptPos0+1)-(min(aThinWinW-1,max(iThin0,10)):0)
		pAccept[mZ,] <- colSums(acceptWindow[curAcceptRows,,drop=FALSE]) / (length(curAcceptRows)*ctrl$thin) 
		temp[mZ,] <- TstepFixed[iThin0*ctrl$thin+ctrl$thin,]
		if( boRecordProposalsIThin ) Y[,iGen+1-nGenOmitRecord,] <- resDo$Y	# recored the proposals together with Density
		
		# update the number of generations to choose from, which depends on acceptance rate
		tmp.chains <- rep(1:nPops, each=nChainsPop)	# population of chain at position
		ar <- pmax(0.026,tapply( pAccept[mZ,], tmp.chains, mean ))	#mean acceptance rate per population, assume minimum 1%=0.16*0.16 to ensure localization (else may go to 0 -> entire history -> no acceptance) 
		nGenBack <- pmin(mZ,ceiling(ctrl$gInd * pmax(1,ctrl$pIndStep/(ar * ctrl$thin))))	#number of genrations to select states from for each population, ctrl$gInd is multiplied by number of rows for one step depending on acceptance rate and thinning but at least one  
		
	} # iThin0 in 0:(nThinnedGen-1)
	#res <- list(parms=Z, rLogDen=rLogDen, thin=ctrl$thin, pAccept=pAccept, call=ccl )
	res <- list(parms=Z, rLogDen=rLogDen, thin=ctrl$thin, pAccept=pAccept, temp=temp, resFLogDenX = chainState$resFLogDenX, Y=Y )
	class(res) <- c( class(res), "twDEMC" )	#monte carlo proposal list
	res
	### list of class \code{twDEMC} (with \code{nStep = M0+nGen%/%thin}) \describe{
	### \item{parms}{ array (d x nStep x nChain) the initial parameters (last row of M0) and the accepted parameter combinations}
	### \item{rLogDen}{ array (nStep x nChain) the logDen of the accepted parameter combinations}
	### \item{pAccept}{ array (nStep x nChain) the acceptance probability over previous ctrl$pAcceptWindowWidth steps }
	### \item{temp}{ vector (nStep) the Temperature used at the step }
	### \item{resFLogDenX}{ array (nInternalResult x nChain) of fLogDen results of currently accepted step}
	### \item{thin}{ numeric: thinning interval }
	### \item{Y}{ numeric matrix: record of all proposals together with results components and columns "rLogDen" and boolean "accepted".
	###		The length is the thinning intervals of the last 128steps or (if doRecordProposals) all the steps +1 for the first line of the initial state.}
	### }
}
#mtrace(twDEMCInt)
attr(twDEMCInt,"ex") <- function(){
	data(twLinreg1); attach(twLinreg1); plot(obs~xval); abline(theta0)

	# setup the model and the function of LogDensity of parameters
	# ?dummyTwDEMCModel
	# ?logDenGaussian
	argsFLogDen <- list( fModel=dummyTwDEMCModel,obs=obs,invCovar=invCovar,	xval=xval )
	do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen)) #test calling the logDen function
	
	# create an initial distribution of states around the estimate of a usual regression 
	lmDummy <- lm( obs ~ xval, weights=1/sdObs^2)		# results without priors
	(.expTheta <- structure(coef(lmDummy),names=c("a","b")) )
	(.expCovTheta <- {tmp<-vcov(lmDummy); dimnames(tmp)<-list(names(.expTheta),names(.expTheta));tmp} )		# a is very weak constrained, negative covariance between a an b
	confint(lmDummy)
	.nPops=2
	Zinit <- initZtwDEMCNormal( .expTheta, .expCovTheta, nChains=4*.nPops, nPops=.nPops)
	
	# run the chains
	res <-  twDEMC( Zinit, nGen=500, fLogDen=logDenGaussian, argsFLogDen=argsFLogDen, nPops=.nPops )
	
	# plot the results
	plotThinned(as.mcmc.list(res))	# no apparent burnin
	
	# get the empirical standard deviation and 95% confidence intervals from the DEMC-sample
	sample <- stackChains(res)
	apply(sample,2, sd)	
	apply(sample,2, quantile, probs=c(0.025,0.975) )
	
	# plot density surface 
	pairs(sample)
	ds <- as.data.frame(apply(sample[ sample[,"rLogDen"] >= max(sample[,"rLogDen"]-1.9), ],2,function(var){
			grain <- diff(range(var))/60
			round(var/grain)*grain
		}))
	tmpx <- sort(unique(ds$a)); tmpy <- sort(unique(ds$b))
	dsog <- expand.grid(x=tmpx[-1], y=tmpy[-1])
	dsog$z <- apply(as.matrix(dsog[,1:2]),1,function(xy){ 
			dss <- subset(ds,a==xy[1] & b==xy[2])
			if( 0<nrow(dss)) max(dss$rLogDen) else NA
		})
	image( tmpx, tmpy,  matrix(dsog$z,nrow=length(tmpx)-1), col = rev(heat.colors(100)), xlab="a", ylab="b" )
	
	#nicer and faster with package lattice
	.tmp.f <- function(){
		# library(lattice)
		# round numbers to see something in levelplot else points get too small
		sampleSig <- apply(sample[ sample[,"rLogDen"] >= max(sample[,"rLogDen"]-1.9), ],2,function(var){
				grain <- diff(range(var))/150
				round(var/grain)*grain
			})
		levelplot(rLogDen~a*b, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	}

	# better seen in three dimensions
	.tmp.f <- function(){
		sampleSig <- sample[ sample[,"rLogDen"] >= max(sample[,"rLogDen"]-1.9), ]
		ds <- as.data.frame(sampleSig)
		plot3d(ds$a, ds$b, ds$rLogDen
		  ,col=rev(heat.colors(100))[round(rescale(ds$rLogDen,to=c(1,100)))]
		  ,xlab="a", ylab="b", zlab="LogDen"
		)
		view3dTiltSpin() 
		play3d( spin3d() )	# press ESC for stop
	}
	
	# for mor examples with prior and simulated annealing see the vignettes
}

.tmp.f <- function(){
	# why
	tmp <- sample[ sample[,"b"]>8 & sample[,"rLogDen"]>-17,,drop=FALSE]
	plot(obs~xval); abline(tmp[1,2:3])
}

.generateXPropThin <- function(
		### Generate Proposal for the next thin steps.
	nPops,
	Z,mZ,
	ctrl,
	nGenBack,
	rLogDen
){
	##seealso<< 
	## \code{\link{twDEMCInt}}
	## \code{\link{.xStepSnooker}}
	## \code{\link{.xStepParallel}}
	
	d <- as.list(structure(dim(Z),names=c("parms","gen","chains")))
	d$gen <- mZ
	d$steps <- ctrl$thin
	nChainsPop = d$chains %/% nPops
	##details<<  
	## Random states for chains for difference vector generation are within subsets of populations.
	## This allows simulating independent population of chains.
	## The acceptance rate may differ amonng populations. Hence, the set of previous generations to 
	## randomly select from also differs between poplations.
	# integer array (thinSteps*nChain*4) sample of chains, within populations
	X <- adrop(Z[,mZ,,drop=FALSE],2)		#assume current state X as beginning of the interval
	nStates <- d$steps*nChainsPop*3
	fzPop <- function(
		### stacked random states and chains in order to vectorize 
		iPop
	){
		sChains <- ((iPop-1)*nChainsPop+1):(iPop*nChainsPop) 
		sGens <- (mZ-nGenBack[iPop]+1):mZ 
		rrGenPop <-  sample(sGens, nStates, replace = TRUE)
		rrChainsPop <-  sample(sChains, nStates, replace = TRUE)
		# in order to constrain two dimensions at one time use the [] subset with an array see ?"["
		rr <- cbind(rrGenPop,rrChainsPop)
		zLogDen <- array(rLogDen[rr], dim=c(1,d$steps*nChainsPop,3), dimnames=list(parms="logDen", steps=NULL, zi=NULL) )
		rrParms <- cbind( rep(1:d$parms, nStates), rep(rrGenPop,each=d$parms), rep(rrChainsPop,each=d$parms) )
		zParms <- array(Z[rrParms], dim=c(d$parms,d$steps*nChainsPop,3), dimnames=list(parms=rownames(Z), steps=NULL, zi=NULL) )
		z0 <- abind( zParms, zLogDen, along=1)
		
		#append as forth initial state x vector to z
		chainsZ <- rep(sChains,each=d$steps)
		Xs <- array(X[,chainsZ], dim=c(d$parms,d$steps*nChainsPop,1), dimnames=list(parms=rownames(Z), steps=NULL, zi=NULL) )
		XsLogDen <- array( rLogDen[mZ,chainsZ],dim=c(1,d$steps*nChainsPop,1), dimnames=list(parms="logDen", steps=NULL, zi=NULL)  ) 
		Xs <- abind( Xs, XsLogDen, along=1)	#logDen row never used for X
		z <- abind( z0, Xs, along=3)
		
		##details<< 
		## States z2 is  attempted to be distinct.
		## And z3 is attempted to be distinct from x (assuming steps in z are stacked chains collectives of X)  
		iSame <- twWhichColsEqual( z[,,2], z[,,1] )
		i <- 0
		while( 0 < length(iSame) && i<5){
			z[,iSame,2] <- z[,sample.int(dim(z)[2],length(iSame)),2]
			iSame <- twWhichColsEqual( z[,,2], z[,,1] )
			i<-i+1
		}
		(tmp <- which( z[,,2] == z[,,1], arr.ind=TRUE))
		#when z3 is the same as x, i.e. z[,,4]
		iSame <- twWhichColsEqual( z[,,3], z[,,4] )
		#iSame <- unique(which( (z[-nrow(z),,3]==Xs), arr.ind=TRUE )[,2])
		i <- 0
		while( 0 < length(iSame) && i<5){
			z[,iSame,3] <- z[,sample.int(dim(z)[2],length(iSame)),3]
			iSame <- twWhichColsEqual( z[,,3], z[,,4] )
			i<-i+1
		}
		(tmp <- which( z[,,3] == z[,,4], arr.ind=TRUE))
		z
		### random states (Nparms+1,(steps*nChainsPop), 4)
		### first dimension is the state vector appended by its logDen
		### random states for each step and chains (stacked to be vectorized)
		### chain is last dimensionion in stack (consequtive steps for one chain) in  order to abind across populations
		### z dim: three random vectors, forth dimension is the initial state x of the chain  
	}
	zxl <- lapply( 1:nPops, fzPop ) 
	zxAndLogDen <- structure( abind(zxl, along=2) , dimnames=list(parms=c(rownames(Z),"logDen"),steps=NULL,zi=NULL))
	zx <- zxAndLogDen[1:d$parms,,,drop=FALSE]
	z <- zx[,,1:3,drop=FALSE]	#three random state vectors per step
	X <- adrop(zx[,,4,drop=FALSE],3)	#initial state vector for step
	zLogDen <- adrop(zxAndLogDen[,,1:3,drop=FALSE][d$parms+1,,,drop=FALSE],1)
	dz <- as.list( structure(dim(z), names=names(dimnames(z))) )
	
	res <- matrix( as.numeric(NA), nrow=dz$parms+1, ncol=dz$steps, dimnames=c(list(c(rownames(z),"rExtra")),list(NULL)) )
	boSnooker <- runif(dz$steps)< ctrl$pSnooker
	if( 0 < sum(boSnooker) ){
		res[,boSnooker] <- .xStepSnooker(z[,boSnooker,,drop=FALSE],X[,boSnooker,drop=FALSE],ctrl=ctrl)
	}
	if( 0 < sum(!boSnooker) )
		res[,!boSnooker] <- .xStepParallel(z[,!boSnooker,,drop=FALSE],zLogDen=zLogDen[!boSnooker,,drop=FALSE],ctrl=ctrl)	
	
	#second dimension is Nsteps*nChains (we can set nChains as last dimension 
	#array(chainsZ,dim=c(d$steps,d$chains))
	resArraySteps <- array(res, dim=c(d$parms+1,d$steps,d$chains), dimnames=list(parms=rownames(res),steps=NULL,chains=NULL) )
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
	
	d <- as.list(structure(dim(z),names=c("parms","steps","iz")))
	gamma_snooker = runif(d$steps, min=1.2,max=2.2)
	res <- matrix( as.numeric(NA), nrow=d$parms+1, ncol=d$steps, dimnames=c(list(c(rownames(z),"rExtra")),list(NULL)) )
	#need loop because of inner products
	for( i in 1:d$steps){
		x_z = X[,i] - z[,i,3]
		D2 = max(1.0e-300, x_z %*% x_z)
		#gamma_snooker =1.7
		projdiff = ((z[,i,1] -z[,i,2]) %*% x_z)/D2  # inner_product of difference with x_z / squared norm x_z
		res[-(d$parms+1),i] <- xStepChain <-  (gamma_snooker[i] * projdiff) * x_z
		xPropChain = X[,i] + xStepChain
		x_z = xPropChain - z[,i,3]
		D2prop = max((x_z%*%x_z), 1.0e-30)
		res[(d$parms+1),i] <- rExtra <- ctrl$Npar12 * (log(D2prop) - log(D2))   # Npar12  =(d$parms - 1)/2  # extra term in logr for accept - reject ratio
	}
	res
} # DE-Snooker update
	
.xStepParallel <- function(
	### DE-parallel direction update based on given random numbers.
	z, ##<< numeric array (Nparms,(nsteps), 3) of random states, dimnames parms,steps,zi
	zLogDen,	##<< numeric matrix (nsteps,3): logDen corresponding to the random states z   
	ctrl
){
	##seealso<< 
	## \code{\link{.generateXPropThin}}
	
	d <- as.list(structure(dim(z),names=c("parms","steps","iz")))
	dz <- adrop((z[,,1,drop=FALSE]-z[,,2,drop=FALSE]),3)	#jump vector as the difference between two random states (from z2 towards z1)
	if( !is.null(ctrl$probUpDir) && (ctrl$probUpDir != 1/2) ){
		iFinites <- which(is.finite(zLogDen[,1]) & is.finite(zLogDen[,2]))
		iUpsF <- which( zLogDen[iFinites,1] < zLogDen[iFinites,2] ) # when proposing a jump between two states towards lower density
		iUps <- iFinites[iUpsF]
		if( 0 < length(iUps)){
			mult <- ifelse( runif(length(iUps)) >= (2*ctrl$probUpDir-1), -1, 1) # increase chance of going directio of upward LogDensity 
			dz[,iUps] <- dz[,iUps] * rep( mult, each=d$parms )
		}
	}
	boGammasF1 <- (runif(d$steps) < ctrl$pGamma1)
	gamma_par <- matrix( ctrl$F1, nrow=d$parms, ncol=d$steps) # to be able to jump between modes
	gamma_par[,!boGammasF1] <- ctrl$F2 * runif( d$parms*sum(!boGammasF1), min=1-ctrl$epsMult, max=1+ctrl$epsMult)    # multiplicative error to be applied to the difference 	
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
	resFLogDenX,	##<< numeric array: result of fLogDen for current location, rows: result components, columns chains  
	logDenX, 		##<< numeric vector of LogDensity of current position of chains
	xStep, 			##<< array rows: difference vectors in parameter space, cols: chains, zdim: steps
	rExtra,			##<< numeric matrix: row: chain, col: step within thinning interval
	temp,			##<< numeric matrix: row: chain, col: step within thinning interval
	nPops,			##<< the number of populations
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
	
	#if(!is.numeric(X) | !is.numeric(resFLogDenX) | !is.numeric(logDenX) )
	#	stop(".doDEMCSteps: first three arguments must be numeric")
	xStepStacked <- do.call( rbind, lapply(1:(dim(xStep)[3]),function(iStep){t(adrop(xStep[,,iStep,drop=FALSE],3))}) )
	#all chains of first step, all chains of second step, ...
	d <- as.list(structure(dim(xStep),names=c("parms","chains","steps"))) 
	iGenT <- (1:d$steps)
	nCases = d$chains * d$steps
	nChainsPop = d$chains / nPops
	if( !(nCases == length(rExtra)) || !(dim(rExtra)==dim(temp)) )
		stop("number of cases in steps must correspond to length of rExtra and length of temp")
	#iPops <- matrix( 1:d$chains, ncol=nPops)	
	#F_APPLY <- .doDEMCStep	#the function executed on the nodes: one metropolis step
	F_APPLY <- sfRemoteWrapper	#the function executed on the nodes: one metropolis step
	#fArgsSame <- list( remoteFunArgs=as.name("argsDEMCStepWrapper") )	#exported, includes remoteFun=.doDEMCStep, remoteDumpfileBasename and argsDEMCStep
	F_ARGS <- function(i,prevRes){ 
			iChain0<-((i-1) %% d$chains)
			iPop<-(iChain0 %/% nChainsPop)+1
			args<-c(	
				prevRes[c("x", "resFLogDenAcc", "logDenAcc")],
				list( step=xStepStacked[i,,drop=TRUE], rExtra=rExtra[i], temp=temp[i], iPop=iPop, pAccept=pAccept),
				#list( argsDEMCStep=argsDEMCStep )
				list( remoteFunArgs=remoteFunArgs )
			)}	
	#.res0 <- lapply(1:nrow(X),function(row){X[row,]})
	res0 <- lapply(1:d$chains,function(iChain){list(
				x=X[,iChain,drop=TRUE]
				,resFLogDenAcc=resFLogDenX[,iChain,drop=TRUE]
				,logDenAcc=logDenX[iChain]
			)})
	res <- sfFArgsApplyDep( nCases, F_ARGS, F_APPLY, res0, debugSequential=debugSequential)
	#all chains of first step, all chains of second step, ...
	#iChains <- (d$steps-1)*d$chains+(1:d$chains)	#index of the chains of the last step
	#modify in place, so that dimnames etc are preserved
	endChain0 <- (d$steps-1)*d$chains	#index before last step of all chains
	#s0 <- list(X=X, resFLogDenX=resFLogDenX, logDenX=logDenX)#save former state
	#boResFLogDenX <- { .nc <- ncol(resFLogDenX); ( !is.null(.nc) && (.nc > 0) ) }	#if number of columns > 0
	#resFLogDenComp
	for( iChain in 1:d$chains ){
		resChain = res[[endChain0+iChain]]
		X[,iChain] <- resChain$x
		resFLogDenX[,iChain] <- resChain$resFLogDenAcc
		logDenX[iChain] <- resChain$logDenAcc
	}
	Y <- NULL
	if( doRecordProposals){ 
		parNames <- rownames(X)
		resCompNames <- names(res[[1]]$resFLogDenProp)
		if( is.null(resCompNames) ) resCompNames <- paste("logDenComp",1:length(res[[1]]$resFLogDenProp),sep="")
		tmp <- c("rLogDen",parNames,"accepted",resCompNames) 
		Y <- array( double(length(tmp)*d$steps*d$chains), dim=c(d=length(tmp),nStep=d$steps,nChain=d$chains),	dimnames=list(comp=tmp,steps=NULL,chains=NULL) )			
		for(iStep in 1:d$steps){
			chain0 <- (iStep-1)*d$chains
			for( iChain in 1:d$chains ){
				i <- chain0+iChain
				Y["accepted",iStep,iChain] <- res[[i]]$accepted
				Y[parNames,iStep,iChain] <- res[[i]]$xProp
				Y[resCompNames,iStep,iChain] <- res[[i]]$resFLogDenProp
				Y["rLogDen",iStep,iChain] <- res[[i]]$logDenProp
			}
		}
	}
	acceptedM <- matrix( sapply(res,function(resChain){ resChain$accepted}), byrow=TRUE, ncol=d$chains )	#rows: steps, cols: chains
	accepted <- colSums(acceptedM)
	
	resDo <- list(	X=X, resFLogDenX=resFLogDenX, logDenX=logDenX, accepted=accepted, Y=Y )
	### list with components \describe{
	### \item{X}{matrix current position, column for each chain}
	### \item{resFLogDenX}{matrix: result components of fLogDen current position, column for each chain}
	### \item{logDenX}{vector current logDen of chains}
	### \item{accepted}{numerical vector, number of accepted steps for each chain}
	### \item{Y}{numerical matrix (steps x 2+params+result): accepted, rLogDen, parms, and all fLogDen result components for each proposal }
	}

.doDEMCStep <- function( 
	### Perfrom one DEMC step, function to be called in remote process.
	x				##<< numeric vector: current state
	,resFLogDenAcc	##<< numeric vector: all components of fLogDen results for current state
	,logDenAcc		##<< numeric scalar: LogDensity of the current state
	#, resFLogDenX	##<< numeric vector: current states components of fLogDen results
	#, logDenXExt	##<< numeric scalar
	,step			##<< numeric vector of step in X space
	,rExtra 		##<< corrector for rLogDen due to selection of step 
	,temp			##<< temperature of the current step and population
	,iPop			##<< the number of the population for which to perform step
	,pAccept		##<< current acceptance rate
	,argsDEMCStep	
### arguments that do not change between steps: list with components \describe{
### \item{fDiscrProp,argsFDiscrProp}{function and additional arguments applied to xProp, e.g. to round it to discrete values}
### \item{argsFLogDen, fLogDenScale}{additional arguments to fLogDen and scalar factor applied to result of fLogDen}
### \item{posLogDenInt}{the matching positions of intResCompNames within the the results components that are handled internally}
### \item{useMultiT}{boolean wheter to scale temperature for different data streams}
### \item{Tprop}{numeric matrix (comp x pops): proportions of temperature for different data streams}
### \item{TFix}{named numeric vector (comp): components with fixed Temperature}
### \item{posTFix}{integer vector (comp): =match(TFix, compNames): positions of TFix within comp provided for performance reasons}
### \item{DRgamma}{ if !0 and >0 delayed Rejection (DR) (Haario06) is applied by jumping only DRgamma distance along the proposal }
### \item{DRgamma}{ if !0 and >0 delayed Rejection (DR) (Haario06) is applied by jumping only DRgamma distance along the proposal }
### }
){
	#.doDEMCStep
	##seealso<< 
	## \code{\link{twDEMCInt}}
	## \code{\link{.doDEMCSteps}}
	#attach( argsDEMCStep )
	#stop(".doDEMCStep: stop to trace error in remote R invocation.")
	with( argsDEMCStep, {
		boResFLogDenX <- (length(posLogDenInt) > 0)
		# LogDensity of accepted state
		La <- resFLogDenAcc*fLogDenScale	#logDensity	of accepted state
		#assume that all is.finite(resFLogDenAcc), make sure in twDEMCInt
		LaExt <- La
		logDenAcc <- sum(La)
		##details<< \describe{\item{Temperature proportions}{
		## If ctrl$useMultiT=TRUE then at high Temperatures, all datastreams are weighted so that 
		## each one has the same influcence.
		## The Temperature of the components with higher LogDensity (less negative)
		## }}
		Ti <- structure(rep(temp,length(La)),names=names(La))
		if( useMultiT & (length(Tprop)>1) ){
			# make sure that all Tprop <= 1
			Ti <- structure( pmax(1, temp*Tprop[,iPop]), names=names(La)) 			
		}
		##details<< \describe{\item{Components with fixed temperature}{
		## If ctrl$TFix is given then the given componenents are assigened the given
		## temperature independent of decreasing global temperature.
		## Useful for priors: \code{TFix = c(parms=1)}
		## }}
		if( 0 < length(posTFix)){			Ti[posTFix] <- TFix		}
		accepted<-FALSE
		xProp = x + step
		if( is.function(fDiscrProp)) xProp = do.call(fDiscrProp,xProp,argsFDiscrProp, quote=TRUE)
		res <- if(boResFLogDenX){
				resFLogDenInt <- resFLogDenAcc[posLogDenInt]
				TiInt <- Ti[posLogDenInt]
				do.call( fLogDen, c(list(xProp, resFLogDenInt, TiInt), argsFLogDen) )	# evaluate logDen
			}else
				do.call( fLogDen, c(list(xProp), argsFLogDen) )	# evaluate logDen
		#take care that the result has always the same sames, even when if fails
		#if( 0==length(names(res)))
		#	stop("encountered result of fLogDen without names")
		#if( !identical(names(resFLogDenAcc),names(res)))
		#	stop("encountered result with different names")
		#strip attributes other than names, else twDynamicClusterApplyDep fails with big data chunks
		attributes(res) <- list(names=names(res))
		logDenProp=-Inf
		Lp <- res*fLogDenScale	# LogDensity of proposal
		#make sure Lp, La have the same order and legnth
		#if( !identical( names(Lp), names(La)) ) stop(".doDEMCStep: resFLogDenAcc must contain the same components and the order of result of fLogDen." )
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
				resFLogDenAcc <- res
			}				
		}else logAlpha10 <- -Inf
		
		if(!accepted && !is.null(argsDEMCStep$DRgamma) && !is.null(argsDEMCStep$minPCompAcceptTempDecr) && (argsDEMCStep$DRgamma > 0) && (pAccept < 1.2*argsDEMCStep$minPCompAcceptTempDecr)) {
			#----- delayed rejection (DR) step
			# only if acceptance rate drops below 1.2*minAcceptrate
			# repeat all above with delayed rejection (DR) step, only adjust DRfac after calculating Lp
			Lp1 <- Lp
			xProp1 <- xProp
			xProp <- x + argsDEMCStep$DRgamma*step
			if( is.function(fDiscrProp)) xProp = do.call(fDiscrProp,xProp,argsFDiscrProp, quote=TRUE)
			res <- if(boResFLogDenX){
					resFLogDenInt <- resFLogDenAcc[posLogDenInt]
					TiInt <- Ti[posLogDenInt]
					do.call( fLogDen, c(list(xProp, resFLogDenInt, TiInt), argsFLogDen) )	# evaluate logDen
				}else
					do.call( fLogDen, c(list(xProp), argsFLogDen) )	# evaluate logDen
			#take care that the result has always the same sames, even when if fails
			#if( 0==length(names(res)))
			#	stop("encountered result of fLogDen without names")
			#if( !identical(names(resFLogDenAcc),names(res)))
			#	stop("encountered result with different names")
			#strip attributes other than names, else twDynamicClusterApplyDep fails with big data chunks
			attributes(res) <- list(names=names(res))
			logDenProp=-Inf
			Lp <- res*fLogDenScale	# LogDensity of proposal
			#make sure Lp, La have the same order and legnth
			#if( !identical( names(Lp), names(La)) ) stop(".doDEMCStep: resFLogDenAcc must contain the same components and the order of result of fLogDen." )
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
					resFLogDenAcc <- res
				}				
			}
		}	## end DR step
		
		#will invoke prevRes[c("x", "resFLogDenAcc", "logDenAcc")]
		if(!is.numeric(x) | !is.numeric(resFLogDenAcc) | !is.numeric(logDenAcc) )
			stop(".doDEMCStep: x, logDenAcc and resFLogDenAcc must be numeric")

		list(accepted=accepted
			,x=x,resFLogDenAcc=resFLogDenAcc, logDenAcc =logDenAcc		# input to repeated call
			,xProp=xProp,resFLogDenProp=res, logDenProp=logDenProp
		)
	})
	#detach( argsDEMCStep ); list(accepted=accepted,x=x,resFLogDenX=resFLogDenX,logDenX=logDenX)
	### list with components \describe{
	### \item{accepted}{boolean scalar: if step was accepted}
	### \item{x}{numeric vector: current position in parameter space}
	### \item{resFLogDenAcc}{numeric vector: result components of fLogDen for current position }
	### \item{logDenAcc}{numeric vector: summed fLogDen for current accepted position}
	### \item{xProp}{numeric vector: proposal}
	### \item{resFLogDenProp}{numeric vector: result components of fLogDen for proposal }
	### \item{logDenProp}{numeric vector: summed fLogDen for proposal}
	### }
}

twCalcLogDenPar <- function(
	### Invokes fLogDen with proposal in a parallel load balanced way.
	fLogDen,				##<< the objective function
	xProp,					##<< numeric matrix of proposals, columns: parameter vector components rows: cases 
	resFLogDenX=NULL	
		### numeric matrix of result of fLogDen
		### colnames must contain intResCompNames 
		### rows: number of cases in xProp	
	,intResCompNames=character(0)	
		### character vector: names of results components of fLogDen that are used for internal Metropolis decisions 
	,argsFLogDen=list()		##<< arguments passed to fLogDen
	,fLogDenScale=1			##<< factor multiplied to the result of fLogDen
	,debugSequential=FALSE	##<< see \code{\link{sfFArgsApplyLB}}
	,remoteDumpfileBasename=NULL,	##<< see \code{\link{sfRemoteWrapper}}
	...						##<< further arguments passed to fLogDen
){
	##seealso<< 
	## \code{\link{twDEMCInt}}
	#if( (0 == length(resFLogDenX)) ) resFLogDenX=xProp[,FALSE,drop=FALSE]
	if( {tmp<-list(...); any(""==names(tmp)) || length(names(tmp))!=length(tmp)} )
		("twCalcLogDenPar: encountered unnamed argument in ... Check for <- and ,, in list()")
	boProvideX2Argument <- (0 < length(intResCompNames))
	if(boProvideX2Argument ){
		if( 0 == length(resFLogDenX) )
			resFLogDenX <- matrix(-Inf*fLogDenScale, ncol=length(intResCompNames), nrow=nrow(xProp), dimnames=list(NULL,parms=intResCompNames))
		if( !is.numeric(resFLogDenX) || !is.matrix(resFLogDenX) || nrow(resFLogDenX)!=nrow(xProp) )
			stop("resFLogDenX must be a numeric matrix with one row for each chain and column names correponding to a subst of names of result vector of fLogDen")
		iNames <- match( intResCompNames, colnames(resFLogDenX) )
		if( any(is.na(iNames)) )
			stop("if resFLogDenX is given, it must contain named columns for each entry of intResCompNames")
		resFLogDenXInt <- resFLogDenX[,iNames,drop=FALSE]		
	}
	res <- if(boProvideX2Argument){
		#call fLogDen with second argument
		F_ARGS <- function(i){c(list(xProp[i,]),list(resFLogDenXInt[i,]))}
		#F_ARGS(1)
		resl <- sfFArgsApplyLB( nrow(xProp), F_ARGS, F_APPLY=sfRemoteWrapper, remoteFun=fLogDen
		, debugSequential=debugSequential, remoteDumpfileBasename=remoteDumpfileBasename, SFFARGSAPPLY_ADDARGS=argsFLogDen, ...) 
		sfSimplifyLBResult(resl)
	}else{
		do.call( sfApplyMatrixLB, c(list( X=xProp, MARGIN=1, FUN=sfRemoteWrapper, remoteFun=fLogDen		, debugSequential=debugSequential, remoteDumpfileBasename=remoteDumpfileBasename), argsFLogDen, list(...)) )	#use doCall in order to use argsFLogDen
	}
	.logDen <- if( is.matrix(res) )
			colSums(res)*fLogDenScale	
		else
			res*fLogDenScale
	.resFLogDen <- if( is.matrix(res) )	t(res)	else matrix(res,ncol=1,dimnames=list(NULL,rownames(res)))
	list( logDen=.logDen, resFLogDen=.resFLogDen)
	### List with the following items \describe{
	### \item{logDen}{numeric vector: for each state: the sum of logDens over all components, multiplied by fLogDenScale}
	### \item{resFLogDen}{numeric matrix: return components of fLogDen, one row for each state, columns: components }
	### }
}
# twUtest("twDEMC","test.twCalcLogDenPar")
attr(twCalcLogDenPar,"ex") <- function(){
	data(twdemcEx1)
	xProp <- stackChains(twdemcEx1$parms)

	data(twLinreg1)
	attach( twLinreg1 )
	res <- twCalcLogDenPar(logDenGaussian,xProp
		,fModel=dummyTwDEMCModel		### the model function, which predicts the output based on theta 
		,obs=obs			### vector of data to compare with
		,invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	str(res)
	resM <- matrix(res$logDen, ncol=dim(twdemcEx1$parms)[3])
	str(resM)
}


setMethodS3("twDEMC","array", function( 
	### Initialize \code{\link{twDEMCInt}} by array of initial population and remove those generations from results afterwards
	Zinit, ##<< initial population: a numeric array (d x M0 x Npop) see details in \code{\link{twDEMCInt}}  
	...	##<< further arguments to \code{\link{twDEMCInt}}
){
	M0 <- dim(Zinit)[2]
	res <- twDEMCInt(Zinit,...)
	#remove the initial M0-1 rows of Zinit
	subset(res, -(1:(M0-1)) )
	### result of \code{\link{twDEMCInt}} with initial M0-1 cases removed
})

setMethodS3("twDEMC","twDEMC", function( 
	### initialize \code{\link{twDEMCInt}} by former run and append results to former run
	Zinit, ##<< list of class twDEMC, result of \code{\link{twDEMCInt}}
	... ##<< further arguments to \code{\link{twDEMCInt}}
){
	.dots <- list(...)
	argsList <- list(Zinit=Zinit$parms)	
	M0 <- nrow(Zinit$rLogDen)
	#extract X, logDenX, and resFLogDenX from Zinit, but only if not given with \dots
	if( is.null(.dots$X) )		#use is.null, because if provided zero length vector, we want to use it 
		argsList$X <- adrop(Zinit$parms[,M0,,drop=FALSE],2)
	if( is.null(.dots$logDenX) ) 
		argsList$logDenX <- Zinit$rLogDen[M0,,drop=TRUE]
	if( is.null(.dots$resFLogDenX) ) 
		argsList$resFLogDenX <- Zinit$resFLogDenX
	if( is.null(.dots$nPops) ) 
		argsList$nPops <- ncol(Zinit$temp)
	res <- do.call( twDEMCInt, c(argsList,.dots))
	res$rLogDen[1:M0,] <- Zinit$rLogDen
	res$pAccept[1:M0,] <- Zinit$pAccept
	res$temp[1:M0,] <- Zinit$temp
	res
	### Zinit appended with the further generations of \code{\link{twDEMCInt}} 
})
#mtrace(twDEMC.twDEMC)
#mtrace(twDEMCInt)
#twUtest(twDEMC,"test.ZinittwDEMC")


#----------------------- DEMCzspBatch
twDEMCBatch <- function(
	### Calls \code{\link{twDEMCBatchInt}} with arguments taken from attribute \code{batchCall} of \code{Zinit}.
	Zinit	##<< the twDEMC object returned by \code{\link{twDEMCInt}}
	,...	
		### Further arguments that are appended/overwrite entries of \code{batchCall}.
		### If Zinit is of class twDEMC, arguemnts  "resFLogDenX","logDenX" are set to NULL, so that these are taken from Zinit itself.
){
	##seealso<< 
	## \code{\link{twDEMCBatchInt}}
	## \code{\link{twDEMCInt}}
	
	if( is.null(attr(Zinit,"batchCall")) )
		##details<< 
		## If twDEMC-list Zinit has no attribute batchCall, then \code{\link{twDEMCBatchInt}} is called with arguments provided.
		res <- twDEMCBatchInt(Zinit,...)
	else{
		##details<< 
		## If twDEMC-list Zinit has attribute \code{batchCall}, then it is re-executed with adjusted Zinit.
		## All further passed arguments will overwrite arguments in batchCall.
		## Hence it is possible to continue a run easily with \code{twDEMCBatch( resultPrev, nGen=nGenPrev+100 )}
		cl <- attr(Zinit,"batchCall")
		if( !is.call(cl)) stop("twDEMCBatch: attr(Zinit,\"batchCall\") is not a call.")
		cl[ c("resFLogDenX","logDenX","X") ] <- NULL	#these arguments will be calculated from Zinit
		#overwrite by further arguments
		cl$Zinit <- Zinit 
		.dots <- list(...)
		cl[ names(.dots) ] <- .dots
		res <- eval(cl)
	}
	res
	### result of \code{\link{twDEMCBatchInt}}
}
#mtrace(twDEMCBatch)

twDEMCBatchInt <- function(
		### Calls \code{\link{twDEMCInt}} successively with \code{nBatch} generations.
	Zinit
		### the initial population
	, nGen=10
		### number of generations in total
	, nBatch=512
		### number of generations between saving results
	, ...
		### further arguments passed to twDEMC
	, restartFilename=NULL
		### name of the file to save intermediate results resRestart.twDEMC, NULL for not saving
	, fCheckConvergence=function(res, addArgs){ FALSE }
		### checking convergence of a DEMC and interrupting 
	, fCheckConvergenceArgs=list()
		### additional arguments to the DEMC convergence checking function
	, T0=1
		### initial temperature of burnin , defaults to 1 or if Zinit is twDEMC to the temperature of the last row
	, nGenBurnin=0
		### initial number of generations of burnin (Temperature decrease to 0, and usage of probUpDirBurnin)
		### either numeric vector of length nPops or scalar
	, doResetOutlierN=256
		### if > 0, outlier chains are reset to best chain
	, nPops = 1
		### number of independent populations
	, probUpDirBurnin=0.75
		### probUbDir during burnin (see twDEMC argument propUpDir)
	, controlTwDEMC=list()
		### controls to \code{twDEMCInt} and twDEMCBatch, some items are overwritten before passing.
		### This depends on convergence and accpetance during burnin.
	#,minPCompAcceptTempDecr=0.2	
		### if maximum di=Lpi-Lai drops below this rate 
		### then Temperature is not decreased within the next batch.
	, doRepeatLowAcceptanceChains=TRUE
	, maxNGenBurnin=50000	##<< maximum burnin beyond which can not be extendend on too low acceptance rate
	, fCalcTStreamDiffLogDen=calcDEMCTempDiffLogDen3	##<< function calculate optimal Temperature and acceptance rates based on Lp-La 
	, argsFCalcTStreamDiffLogDen=list()		##<< further arguments to fCalcTStreamDiffLogDen  
	, fCalcTGlobal=calcDEMCTempGlobal2		##<< function calculating global target temperature
	, argsFCalcTGlobal=list()				##<< further arguments to fCalcTStreamDiffLogDen  
){
	# twDEMCBatchInt
	##seealso<<   
	## \code{\link{twDEMCInt}}
	## \code{\link{twDEMCBatch}}
	##details<< 
	## Usually invoked by \code{\link{twDEMCBatch}}.
	if( !hasArg(nBatch))nBatch=512
	if( !hasArg(restartFilename))restartFilename=NULL
	cl <- match.call()
	for( i in (2:length(cl)) )
		try( cl[i] <- list(eval.parent(cl[[i]])) )
	#cl[ names(cl)[-1] ] <- lapply(as.list(cl)[-1],eval.parent) 	#substitute all variables by their values, -1 needed for not substituting the function name
	iRun <- if( is(Zinit,"twDEMC") ) getNGen(Zinit) else 0   #already completed generations
	if( length(nGenBurnin)==1 ) nGenBurnin=rep(nGenBurnin,nPops)    # maybe different for each generation
	newNGenBurnin <- nGenBurnin	#initialize
	nRun <- min(nBatch, (if(iRun<max(nGenBurnin)) min(max(nGenBurnin),nGen) else nGen) -iRun)
	
	##details<< \describe{ \item{Detailed control parameters: \code{controlTwDEMC}}{
	##describe<<  
	ctrl = list( 
		minPCompAcceptTempDecr=0.2  ##<< if acceptance rate drops below minPCompAcceptTempDecr+0.02 times this level, employ delayed rejection (DR)
		,nSampleGelmanDiag=512	
			### sample size of the appropriately thinned sample that will be detrended to check Gelman diagnostics before decreasing temperature.
			### When working in high dimensional parameter space, you may encounter false high gelman diagnostics. Then you can try increasing this parameter,
			### but be aware of increasing computational demand.
	)  ##end<< 
	## }}
	ctrl[names(controlTwDEMC)] <- controlTwDEMC
	
	ctrl <- controlTwDEMC
	ctrl$T0=T0
	ctrl$Tend=T0 	#no Temp decrease in first batch (different for Zinit is not twDEMC, see below) 
	ctrl$probUpDir=(if(nRun<=max(nGenBurnin)) probUpDirBurnin else NULL)
	TFix <- if(is.numeric(ctrl$TFix)) ctrl$TFix else numeric(0) 
	thin <- if(is.numeric(ctrl$thin)) 	ctrl$thin else 1
	nGen <- (nGen %/% thin)*thin

	
	##details<< 
	## If Zinit is of class twDEMC, initial temperature is set to the temperature of the last row
	## and the number of generations already in Zinit are skipped.
	pTarget=ctrl$minPCompAcceptTempDecr+0.02
	if( is(Zinit,"twDEMC") ){
		res <- Zinit
		iRun <- getNGen(res)
		if( iRun >= nGen ) return(res)
		ctrl$initialAcceptanceRate <- popMeansTwDEMC( res$pAccept[nrow(res$pAccept),],nPops )
		ctrl$T0 <- T0c <- res$temp[ nrow(res$temp), ,drop=FALSE ]
		if( length(T0c) != nPops) stop(paste("twDEMCInt: encoutered temperature recored with",length(T0c),"columns but argument nPops=",nPops))
		
		#calculate optimal end temperature
		resCols <- match( .getResFLogDenNames(res$resFLogDen), rownames(res$Y))
		#nPops <- ncol(res$temp)
		nChains <- dim(res$parms)[3]
		nChainsPop <- nChains %/% nPops
		chain2Pop <- rep(1:nPops, each=nChainsPop )	#mapping of chain to population
		
		#diffLogDen <- getDiffLogDen.twDEMCProps(res$Y, resCols, nLastSteps=ceiling(128/nChainsPop)) 	#in twDEMC S3twDEMC.R
		diffLogDen <- getDiffLogDen.twDEMCProps(res$Y, resCols, nLastSteps=128) 	#in twDEMC S3twDEMC.R
		diffLogDenPops <- popApplyTwDEMC( diffLogDen, nPops=nPops, function(x){ abind(twListArrDim(x),along=2,new.names=dimnames(x)) })	#stack param columns by population
		Ti <- matrix(1,nrow=nrow(res$resFLogDenX),ncol=nPops, dimnames=list(comp=.getResFLogDenNames(res$resFLogDenX),pop=NULL)) #temperature by component x population 
		pAcceptTVar <- numeric(nPops)
		#newNGenBurnin <- res$nGenBurnin #keep previous, do not reinitialize: res from twDEMC has no nGenBurnin
		TGlobal <- T0c
		for( iPop in which( TGlobal>1 )){
			resPop <- subChains(res,iPops=iPop) 
			dLp=adrop(diffLogDenPops[,,iPop,drop=FALSE],3)
			Ti[,iPop] <- TiPop <- do.call( fCalcTStreamDiffLogDen, c(list(diffLogDen=dLp,TFix=TFix,Tmax=T0c[iPop],pTarget=pTarget), argsFCalcTStreamDiffLogDen) ) # optimal Temperature estimated by dLp for each variable
			maxTiPop <- max(TiPop)	
			pAcceptTVar[iPop] <- attr(TiPop,"pAcceptTVar") # acceptance rate of the Temperature dependent step
			tmp <- do.call( fCalcTGlobal, c(list(resPop=resPop,diffLogDen=dLp,TLp=maxTiPop,pAcceptTVar=pAcceptTVar[iPop],iRun=iRun, nGenBurnin=nGenBurnin[iPop], nRun=nRun),argsFCalcTGlobal) )
			TGlobal[iPop] <- tmp$TGlobal
			newNGenBurnin[iPop] <- tmp$nGenBurnin
		}
		ctrl$Tprop <- Ti		# will be scaled in twDEMC
		ctrl$Tend <- TGlobal			# standard is using Multi-temperature
		nGenBurnin <- pmin(maxNGenBurnin, newNGenBurnin)
		cat(paste("  burnin=",paste(nGenBurnin,collapse="  "),"\n",sep=""))			
		#recalculate nRun with possibly changed nGenBurnin
		nRun <- min(nBatch, (if(iRun<max(nGenBurnin)) min(max(nGenBurnin),nGen) else nGen) -iRun)
		
		##details<< \describe{\item{Temperature estimate from proposal distribution}{
		## The distribution of differences between Density of proposals Lp and of accepted state La
		## can be used to estimate an optimal temperature per data stream, so that each
		## datastream contributes to rejections in about the same magnitude and the overall
		## acceptance rate is about a specified value.
		## The proportions of the so calculated datastream specific temperature are multiplied 
		## with the global temperature on Metropolis decisions.
		## deprecated[Further, if the temperatures of the datasteams are all below the goal of the global
		## temperature Tend, Tend is also lowered]
		## }}
	}
	
	cat(paste(iRun," out of ",nGen," generations completed. T=",paste({T<-ctrl$T0;round(T,digits=ifelse(T>20,0,1))},collapse="  "),"     ",date(),"\n",sep=""))
	
	#--------- do the twDEMC ------
	nRun <- (nRun%/%thin)*thin  #make nRun multiple of thin
	res <- twDEMC( Zinit=Zinit, nGen=nRun, nPops=nPops, controlTwDEMC=ctrl, ... )
	res$nGenBurnin <- newNGenBurnin	
	attr(res,"batchCall") <- cl
	boConverged=FALSE
	if( hasArg(fCheckConvergence))
		boConverged = (all(res$temp[nrow(res$temp),]<1.1)) & fCheckConvergence(res, fCheckConvergenceArgs)
	
	#iRun <- iRun + nRun		#current number of runs # because of thinning may actually performed fewer runs than nRun
	iRun <- getNGen(res)
	nChains <- dim(res$parms)[3]
	nChainsPop <- nChains %/% nPops
	chain2Pop <- rep(1:nPops, each=nChainsPop )	#mapping of chain to population
	resCols <- match( twDEMC:::.getResFLogDenNames(res$resFLogDen), rownames(res$Y))	#index of columns of results components in Y
	while( !boConverged & (iRun < nGen) ){
		cat(paste(iRun," out of ",nGen," generations completed. T=",paste({T<-res$temp[nrow(res$temp),];round(T,digits=ifelse(T>20,0,1))},collapse="  "),"     ",date(),"\n",sep=""))
		##details<< \describe{\item{Saving and restarting}{ 
		## If \code{restartFilename} has been specified then result are stored as variable \code{resRestart.twDEMC} to the given file.
		## Runs can be continued without the need of respecifying all the parameters (see example). 
		## This is helpful with expensive models and long cluster runs for cases where the program has to be aborted.
		## }}
		if(is.character(restartFilename)){
			resRestart.twDEMC = res #avoid variable confusion on load by giving a longer name
			resRestart.twDEMC$nGenBurnin <- newNGenBurnin	# also store updated calculation of burnin time
			save(resRestart.twDEMC, file=restartFilename)
			cat(paste("Saved resRestart.twDEMC to ",restartFilename,"\n",sep=""))
		}
		nRunPrev <- nRun
		zGen <- dim(res$parms)[2]
		# do not remove outliers after burnin
		if( (doResetOutlierN>0) & (iRun <= max(nGenBurnin)) ){
			iGenOmega <- max(1,zGen-doResetOutlierN+1):zGen #(zGen%/%2):zGen
			# according to Vrugt09
			omega <- sapply( 1:nChains, function(iChain){mean(res$rLogDen[iGenOmega,iChain], na.rm=TRUE)}) #mean logDen across last halv of chain
			for( iPop in 1:nPops ){
				cat(paste("     pop ",iPop,": ",sep=""))
				iChains <- ((iPop-1)*nChainsPop+1):(iPop*nChainsPop)
				q13 <- quantile( omega[iChains], c(1/4,3/4) )	#lower and upper quartile
				bo <- omega[iChains] < q13[1] -2*diff(q13)		#outside 2 interquartile ranges, see Vrugt09
				if( any(bo) ){
					cat(paste("   resetting outlier chains ",paste(iChains[bo],collapse=","),"of population ",iPop,"\n"))
					#reset state of outliers to the best sampled parameter state
					tmp.best <- which( res$rLogDen[,iChains] == max(res$rLogDen[,iChains]), arr.ind = TRUE )[1,]	
					res$parms[,zGen,iChains[bo]] <- res$parms[,tmp.best[1],iChains[ tmp.best[2] ] ]
				}
			}
		}
		nRun <- min(nBatch, (if(iRun<max(nGenBurnin)) min(max(nGenBurnin),nGen) else nGen) -iRun)		#iRun: Generation after batch run
		.dots <- list(...)
		.dots[c("logDenX","resFLogDenX")] <- NULL;	#those will be inferred from res
		clArgs <- c(list(Zinit=res), .dots)	#Zinit must be first argument 
		clArgs$nPops<-nPops
		#clArgs$T0=max(1,b*exp(-a*iRun))		# if temp did not decrease start from this temperature
		clArgs$controlTwDEMC <- controlTwDEMC
		clArgs$controlTwDEMC$initialAcceptanceRate <- popMeansTwDEMC( res$pAccept[nrow(res$pAccept),],nPops )
		clArgs$controlTwDEMC$T0<-T0c<-res$temp[ nrow(res$temp), ,drop=FALSE]
		#clArgs$Tend=max(1,b*exp(-a*(iRun+nRun)))
		##--calculating end temperature
		clArgs$controlTwDEMC$Tend <- 1
		if((iRun+nRun)<max(nGenBurnin)) {
			diffLogDen <- getDiffLogDen.twDEMCProps(res$Y, resCols, nLastSteps=ceiling(128/nChainsPop)) 	#in twDEMC S3twDEMC.R
			diffLogDenPops <- popApplyTwDEMC( diffLogDen, nPops=nPops, function(x){ abind(twListArrDim(x),along=2,new.names=dimnames(x)) })	#stack param columns by population
			
			##details<< \describe{\item{cooling and acceptance rate}{ 
			## If acceptance rate of some population drops below rate=ctrl$minPCompAcceptTempDecr then cooling is too fast.
			## In this moderate case do not repeat the rund but keep the current temperature for the next period for this population.
			## and extend the burnin phase by the length of this period.
			## }}
			Ti <- matrix(1,nrow=nrow(res$resFLogDenX),ncol=nPops, dimnames=list(comp=twDEMC:::.getResFLogDenNames(res$resFLogDenX),pop=NULL))
			pAcceptTVar <- numeric(nPops)
			#newNGenBurnin <- res$nGenBurnin # use old one, do not re-initialize, res from twDEMCInt has no entry nGenBurnin
			TGlobal <- T0c
			for( iPop in which( TGlobal>1 )){
				cat(paste("     pop ",iPop,": ",sep=""))
				resPop <- subChains(res,iPops=iPop) 
				dLp=adrop(diffLogDenPops[,,iPop,drop=FALSE],3)
				Ti[,iPop] <- TiPop <- do.call( fCalcTStreamDiffLogDen, c(list(diffLogDen=dLp,TFix=TFix,Tmax=T0c[iPop],pTarget=pTarget), argsFCalcTStreamDiffLogDen) ) # optimal Temperature estimated by dLp for each variable
				maxTiPop <- max(TiPop)	
				pAcceptTVar[iPop] <- attr(TiPop,"pAcceptTVar") # acceptance rate of the Temperature dependent step
				tmp <- do.call( fCalcTGlobal, c(list(resPop=resPop,diffLogDen=dLp,TLp=maxTiPop,pAcceptTVar=pAcceptTVar[iPop],iRun=iRun, nGenBurnin=nGenBurnin[iPop], nRun=nRun),argsFCalcTGlobal) )
				TGlobal[iPop] <- tmp$TGlobal
				newNGenBurnin[iPop] <- tmp$nGenBurnin
			}
			clArgs$controlTwDEMC$Tprop <- Ti		# will be scaled in twDEMC
			clArgs$controlTwDEMC$Tend <- TGlobal			
			nGenBurnin <- pmin(maxNGenBurnin, newNGenBurnin)
			#cat(paste("burnin=",paste(nGenBurnin,collapse="  "),"T=",paste(round(Ti),collapse="  "),"\n",sep=""))			
			cat(paste("  burnin=",paste(nGenBurnin,collapse="  "),"\n",sep=""))			
			#recalculate nRun with possibly changed nGenBurnin
			nRun <- min(nBatch, (if(iRun<max(nGenBurnin)) min(max(nGenBurnin),nGen) else nGen) -iRun)
		}
		nRun <- (nRun%/%thin)*thin  #make nRun multiple of thin
		clArgs$nGen<-nRun
		clArgs$controlTwDEMC$probUpDir <- (if((iRun+nRun)<=max(nGenBurnin)) probUpDirBurnin else NULL)	#set to NULL after burnin
		res <- do.call( twDEMC, clArgs, quote=TRUE )
		res$nGenBurnin <- newNGenBurnin
		attr(res,"batchCall") <- cl
		#res <- twDEMC( Zinit=res, nGen=nRun, ... ) problems with double Zinit 
		if( hasArg(fCheckConvergence))
			boConverged = (all(res$temp[nrow(res$temp),]<1.1)) & fCheckConvergence(res, fCheckConvergenceArgs)
		iRun <- getNGen(res)
	}
	cat(paste(iRun," out of ",nGen," generations completed. T=",paste({T<-res$temp[nrow(res$temp),];round(T,digits=ifelse(T>20,0,1))},collapse="  "),"     ",date(),"\n",sep=""))
	res$nGenBurnin <- newNGenBurnin
	res
	### List of class twDEMC (see \code{\link{twDEMCInt}}) with additional entry nGenBurnin
}
#mtrace(twDEMC)
#mtrace(twDEMCBatch)
#mtrace(twDEMCBatchInt)
attr(twDEMCBatchInt,"ex") <- function(){
	data(twLinreg1); attach( twLinreg1 ) 
	
	# run with saving a restart file each 50 generations
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,				### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior= thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=4*.nPops, nPops=.nPops)
	restartFilename="exampleTwDEMCBatch_saveAndRestart.RData"
	unlink(restartFilename)
	res <-  twDEMCBatch( Zinit=Zinit, nGen=60,	nPops=.nPops
		, fLogDen=logDenGaussian, argsFLogDen=argsFLogDen
		, nBatch=50, restartFilename=restartFilename	# save each 50 generations
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
























logDenGaussian <- function(
	### Invokes the model and calculate an unnormalized logDensity (-1/2*misfit) assuming (multivariate) Gaussian errors in both data in priors. 
	theta,			##<< the parameter vector.
	logDenAccept=numeric(0),	##<< scalar: logDen for parms from revious run for two step Metropolis decision
	metropolisStepTemp=c(parms=1),		##<< numeric named vector: the temperature for internal metropolis step
	..., 			##<< any other arguments passed to fModel
	fModel,			##<< the model function, which predicts the output based on theta 
	theta0=theta,	##<< parameter vector, first argument to fModel. Before invocation components theta overwrite theta0 
	obs,			##<< vector of data to compare with
	invCovar,		##<< the inverse of the Covariance of obs (its uncertainty)
	thetaPrior = NULL,	##<< the prior estimate of the parameters
	invCovarTheta = NULL,	##<< the inverse of the Covariance of the prior parameter estimates
	namesTheta=NULL, ##<< names assigned to theta (if not NULL), before invoking mofModel	
	scale=-1/2 	 		##<< factor to mulitply the misfit (e.g. -1/2 to obtain the unnormalized logDensity)
){
	# logDenGaussian
	##seealso<<   
	## \code{\link{twDEMCInt}}
	## \code{\link{dummyTwDEMCModel}}
	if( !is.null(namesTheta ))
		names(theta) <- namesTheta	##details<<
	theta0[names(theta)] <- theta
	##details<< 
	## If thetaPrior is not specified (NULL) then no penalty is assigned to parameters.
	logDenPropParms <- if( !is.null(thetaPrior) ){
		tmp.diffParms <- theta0 - thetaPrior
		as.numeric(t(tmp.diffParms) %*% invCovarTheta %*% tmp.diffParms) 
	} else 0
	##details<<
	## Supports a two-step Metropolis descision. If \code{logDenAccept["parms"]} is provided, 
	## then a Metropolis descision is done based only on the parameters.
	## If it fails, then \code{c(obs=NA, parms=-Inf)} is returned. 
	## The possible costly evaluation of fModel is avoided.
	if( !is.na((logDenXParms <- logDenAccept["parms"])) & (logDenXParms<0)){
		logr = (scale*logDenPropParms - logDenXParms) / metropolisStepTemp["parms"]
		if ( is.numeric(logr) & (logr) <= log(runif(1)) ){
			#reject
			return(c(obs=NA, parms=-Inf))
		}
	}
	# evaluate the model at parameters theta0 
	tmp.pred <- fModel(theta0, ...)
	tmp.diffObs <- tmp.pred - obs
	tmp.misfit <-  c(obs=as.numeric(t(tmp.diffObs) %*% invCovar %*% tmp.diffObs), parms=logDenPropParms)
	scale * tmp.misfit
	### the misfit: scale *( t(tmp.diffObs) %*% invCovar %*% tmp.diffObs + t(tmp.diffParms) %*% invCovarTheta %*% tmp.diffParms )
}
#mtrace(logDenGaussian)
#mtrace.off()
#twUtestF(logDenGaussian)

dummyTwDEMCModel <- function(
		### example model function: y=a+bx
	theta,	##<< parameter vector with names a and b
	xval	##<< additional argument, passed by ... in logDenGaussian
){ 
	# dummyTwDEMCModel
	##seealso<<   
	## \code{\link{twDEMCInt}}
	## \code{\link{logDenGaussian}}
	theta["a"] + theta["b"]*xval 
}

twRunDEMC <- function(
	### Run a twDEMC with first replacing nonfinite Densitys in Zinits last row.
	#Zinit				##<< initial values
	#,nPops				##<< number of populations
	#,fLogDen			##<< objective function
	#,argsFLogDen		##<< further arguments to objective function
	#,resFLogDenX=NULL	## character vector: components of result of fLogDen that are handled with internal Metropolis steps.
	argsTwDEMCBatch	 
		### Arguments passed to twDEMCBatch -> twDEMCInt -> fLogDen.
		### It is updated by \dots.
		### After update it must contain entries Zinit and fLogDen
		### It is further searched for entries nPops, and argsFLogDen, and resFLogDenX. The latter are initialized to defaults  \code{1,list(),character(0)} respectively if not found.   
	,...				 ##<< further arguments passed to twDEMCBatch -> twDEMCInt -> fLogDen
	,argsReplaceZinit=list() ##<< further arguments passed to twDEMCBatch
	,prevResRunCluster=NULL	##<< results of call to twRunDEMC, argument required to be called from runCluster.R
	,restartFilename=NULL	##<< name of the file to store restart information, argument required to be called from runCluster.R 
){
	#update argsDEMC to args given 
	argsDEMC <- argsTwDEMCBatch
	
	.dots <- list(...)
	argsDEMC[ names(.dots) ] <- .dots
	#for( argName in names(.dots) ) argsDEMC[argName] <- .dots[argName]
	if( !is.null(restartFilename)) argsDEMC$restartFilename <- restartFilename
	
	# get entries from argsDEMC
	Zinit <- argsDEMC$Zinit
	if( is.null(Zinit) ) stop("must provide Zinit with ... or with argsTwDEMCBatch")
	argsDEMC$Zinit <- NULL	# Zinit put separately as first argument in do.call
	if( is(prevResRunCluster,"twDEMC") ){
		#mtrace(twDEMCBatch)
		do.call( twDEMCBatch, c( list(Zinit=prevResRunCluster), argsDEMC) )  	
	}else if( is(Zinit,"twDEMC") ){
		do.call( twDEMCBatch, c( list(Zinit=Zinit), argsDEMC) )  	
	}else{
		fLogDen <- argsDEMC$fLogDen
		if( is.null(fLogDen) ) stop("must provide fLogDen with ... or with argsTwDEMCBatch")
		argsFLogDen <- argsDEMC$argsFLogDen
		if( is.null(argsFLogDen) ) argsFLogDen=list()
		resFLogDenX <- argsDEMC$resFLogDenX
		if( (0 == length(resFLogDenX)) ) resFLogDenX<-character(0)
		nPops <- argsDEMC$nPops
		if( is.null(nPops) ) nPops=1
		
		if( 0<length(argsDEMC$debugSequential) && argsDEMC$debugSequential ) 
			argsReplaceZinit$debugSequential=TRUE
		if( 0<length(argsDEMC$stopOnError) && argsDEMC$stopOnError ) 
			argsReplaceZinit$stopOnError=TRUE
		ZinitLogDen <- do.call( replaceZinitNonFiniteLogDensLastStep, c( list(Zinit,fLogDen=fLogDen,nPops=nPops
			, argsFLogDen=argsFLogDen, resFLogDenX=resFLogDenX), argsReplaceZinit) )
		argsDEMC <- within( argsDEMC, {
				X <- adrop(ZinitLogDen$Zinit[,ncol(ZinitLogDen$Zinit),,drop=FALSE],2)
				logDenX<-ZinitLogDen$rLogDen
				resFLogDenX<-if(0==length(ZinitLogDen$resFLogDen)) character(0) else t(ZinitLogDen$resFLogDen)
			})
		res <- do.call( twDEMCBatch, c(list(Zinit=ZinitLogDen$Zinit),argsDEMC) )
		res
		### results of \code{\link{twDEMCBatch}}
	}
}

twRunFLogDenPar <- function(
	### Wrapper for \code{\link{twCalcLogDenPar}} to accept arguments provided by runCluster.R
	...				 		##<< further arguments passed \code{\link{twCalcLogDenPar}}
	,prevResRunCluster=NULL	##<< results of call to twRunDEMC, argument required to be called from runCluster.R
	,restartFilename=NULL	##<< name of the file to store restart information, argument required to be called from runCluster.R 
){
	twCalcLogDenPar(...)
}

.getResFLogDenNames <- function(
	### Extracting/Creating names for result components of resFLogDen
	resFLogDenX
){
	if( is.array(resFLogDenX)){
		res <- rownames(resFLogDenX)
		return( if( !is.null(res) ) res else paste("logDenResComp",1:nrow(resFLogDenX),sep="") )
	}else{
		res <- names(resFLogDenX)
		return( if( !is.null(res) ) res else paste("logDenResComp",1:length(resFLogDenX),sep="") )
	}
}








