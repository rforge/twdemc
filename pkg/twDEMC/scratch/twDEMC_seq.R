# code from twDEMC before vectorization

.generateXPropThinSeq <- function(
	nPops,
	Z,
	ctrl,
	...
){
	d <- as.list(structure(dim(Z),names=c("parms","steps","chains")))
	d$steps <- ctrl$thin
	nChainPop = d$chains %/% nPops
	# integer array (thinSteps*nChain*3) sample of chains, within populations
	rrcPar <- abind( lapply( 1:nPops, function(i){
				tmp <- ((i-1)*nChainPop+1):(i*nChainPop)
				rrcPar <- array( sample(tmp, d$steps*nChainPop*3, replace = TRUE), dim=list(gen=d$steps,chain=nChainPop,i=3) )
			}), along=2 )
	xStepAndExtraL <- lapply(1:d$steps, function(iGenThin){ 
			#.generateXPropChains(Z=Z,nChains=nChains,X=X,rLogDen=rLogDen,mZ=mZ,Npar=Npar,iGen=iGen,ctrl=ctrl,nChainPop=nChainPop,g=g,rrcPar=rrcPar) 
			.generateXPropChains(iGenThin, Z=Z,ctrl=ctrl,rrcPar=rrcPar,nChainPop=nChainPop,d=d,...) 
		}) #rows: difference vectors in parameter space, cols: chains
	xStepAndExtra <- abind( xStepAndExtraL, rev.along=0)		#third dimension is the step within thinning interval
	xStep <- xStepAndExtra[-nrow(xStepAndExtra),,,drop=FALSE]			#dito
	rExtra <- adrop(xStepAndExtra[nrow(xStepAndExtra),,,drop=FALSE],1)			#second dim (columns step within Thinning interval)
	list(xStep=xStep, rExtra=rExtra)
	### a list with components (nSeps=thin) \describe{
	### \item{xStep}{numeric array (Npar,Nchain,Nsteps): difference vectors in parameter space}
	### \item{rExtra}{numeric matrix (Npar,Nsteps): some extra LogDensity from snooker update}}
}

.xStepParallel.seq <- function(
	### non-vectorized version of xStepParallel for testing proper vectorization
	z, ##<< numeric array (Nparms,(nsteps), 3) of random states, dimnames parms,steps,zi
	#X, ##<< current state (Nparms,(nsteps)) corresponding to chain of second dimension in z
	zLogDen,	##<< numeric matrix (nsteps,3): logDen corresponding to the random states z   
	ctrl
){
	d <- as.list(structure(dim(z),names=c("parms","steps","iz")))
	zOrig <- z	#original z with second step dimension
	
	#generate random numbers in the same sequence as in vectorized version
	boGammasF1 <- (runif(d$steps) < ctrl$pGamma1)
	gamma_par_unif <- runif( d$parms*sum(!boGammasF1), min=1-ctrl$epsMult, max=1+ctrl$epsMult)    # multiplicative error to be applied to the difference
	epsAdd_rnorm <- rnorm(length(dz),0,ctrl$epsAdd)
	
	fStep <- function(iStep){
		z <- zOrig[,iStep]
		dz <- (z[,1]-z[,2])	#jump vector as the difference between two random states (from z2 towards z1)
		if( !is.null(ctrl$probUpDir) && (ctrl$probUpDir != 1/2) ){
			lz <- sapply( 1:2, function(zi){ rLogDen[ rrGen[zi], rrcPar[iGenThin,iChain,zi] ] })#get the logDensity for the random states
			if( all(is.finite(lz)) && (lz[1] < lz[2]) )	# when proposing a jump between two states towards lower density
				if( runif(1) >= (2*ctrl$probUpDir-1) ) dz = -dz	# increase chance of going directio of upward LogDensity
		}
		if ( runif(1)< ctrl$pGamma1 ) { 
			gamma_par = ctrl$F1 # to be able to jump between modes
		} else {
			gamma_par = ctrl$F2 * runif(d$parms, min=1-ctrl$epsMult, max=1+ctrl$epsMult)    # multiplicative error to be applied to the difference
			# gamma_par = F2 
		}
		if (ctrl$epsAdd ==0) {  # avoid generating normal random variates if possible
			#xPropChain = X[,iChain] + gamma_par * dz 
			xStepChain = gamma_par * dz 
		} else {
			#xPropChain = X[,iChain] + gamma_par * dz  + rnorm(d$parms,0,ctrl$epsAdd)
			xStepChain = gamma_par * dz  + rnorm(d$parms,0,ctrl$epsAdd)
		}
		xStepChain
	}
	resl<- lapply( 1:d$steps, fStep)
	xStep <- abind( resl, 2)
	xStepChainAndRExtra <- rbind(xStepChain,rExtra=0)
	### will not give the same result as the parallel version, because rnorm and runif are called for each step
}

.generateXPropChains <- function(
	### generates a proposal for next step for all chains
	iGenThin,	##<< the generation within thinning interval	
	Z,			##<< 
	d,			##<< list of dimensions of the array of steps to generate
	...			##<< arguments passed to \code{\link{.generateXProp}}
){
	xStepAndExtra <- matrix(NA, ncol=dim(Z)[3], nrow=dim(Z)[1]+1, dimnames=list( parms=c(rownames(Z),"rExtra") ,chains=NULL) )
	for (iChain in 1:d$chains){ #for each chain
		xStepAndExtra[,iChain] <- .generateXProp(iGenThin, iChain,Z=Z,d=d,...)
	}#iChain Chains gnerate proposals
	xStepAndExtra
	### numeric matrix
	### each column corresponds to a chain
	### first components of the vector correspond to proposal and last component to some extra logDensity
}

.generateXProp <- function(
	### generates a proposal for next step for chain iChain
	iGenThin,	##<< the generation within thinning interval
	iChain,		##<< the chain number
	Z,			##<< numeric array, space for all parameters for all generations
	rLogDen,	##<< numeric matrix, space for calculated LogDensity
	mZ,			##<< integer, completed generations im mZ and rLogDen
	#Npar,		##<< number components of parameter vector (passed for performace reason instead of calling length or dim)
	#X,			##<< numeric vector, the currently accepted parameter set
	#iGen,		##<< the current generation
	ctrl,  		##<< list twDEMC control parameters
	d,			##<< list of dimensions of step array to create
	nChainPop,	##<< the number of populations
	g,			##<< the number of generations to select past states from
	rrcPar		##<< integer array (nGen*nChain*3) sample of chains, sampled in one step before for performance reasons
){
	X <- Z[,mZ,]		#replace current state X by state of the beginning of thinning interval
	iPop0 = (iChain-1)%/%nChainPop; iPop=iPop0+1 #+1 counting from zero
	rrGen <- sample( (mZ-g[iPop]+1):mZ, 3, replace=FALSE )	#depends on acceptance ratio of population, so cannot optimize performance by sampling togehter beforehand
	if ( runif(1)< ctrl$pSnooker ) {  
		#zChains <- Z[,,rrcSnooker[iGenThin,iChain,] ]	# select the chains
		rrChains <- (iPop0)*nChainPop + sample( nChainPop, 3, replace = TRUE ) #for snooker update all chains of population are equal
		z <- sapply( 1:3, function(zi){Z[,rrGen[zi],rrChains[zi] ]})
		#make sure z[,1] and z[,2] are different states (iChain.e chain did not stay at one place), and also that X[,iChain] and z[,3]] are different
		tmp.i <- 0
		while( (all(z[,1] == z[,2]) | all(X[,iChain] == z[,3])) & (tmp.i < 20)  ){
			rrGen <- sample( (mZ-g[iPop]+1):mZ, 3, replace=FALSE )	
			z <- sapply( 1:3, function(zi){Z[,rrGen[zi], rrChains[zi] ]})
			tmp.i = tmp.i + 1
		}
		# DE-Snooker update
		x_z = X[,iChain] - z[,3]
		D2 = max(sum(x_z*x_z),1.0e-300)
		gamma_snooker = runif(1, min=1.2,max=2.2)
		#gamma_snooker =1.7
		projdiff = sum((z[,1] -z[,2]) *x_z)/D2  # inner_product of difference with x_z / squared norm x_z
		xStepChain = (gamma_snooker * projdiff) * x_z
		xPropChain = X[,iChain] + xStepChain
		x_z = xPropChain - z[,3]
		D2prop = max(sum(x_z*x_z), 1.0e-30)
		rExtra = ctrl$Npar12 * (log(D2prop) - log(D2))   # Npar12  =(d$parms - 1)/2  # extra term in logr for accept - reject ratio
	} # DE-Snooker update
	else {
		# DE-parallel direction update
		# select to random different individuals from Z
		# higher probability to sample from curent chain and that second sample is from same chain as first one 
		#zChains <- Z[,,rrcPar[iGenThin,iChain,] ]	# select the chains
		#z <- sapply( 1:2, function(zi){zChains[,rrGen[zi],zi]})
		z <- sapply( 1:2, function(zi){Z[,rrGen[zi],rrcPar[iGenThin,iChain,zi]]})
		#make sure z[,1] and z[,2] are different states (i.e chain did not stay at one place), and also X[,iChain] and z[,3]]
		tmp.i <- 0
		while( all(z[,1] == z[,2]) & (tmp.i < 20) ){
			rrGen <- sample( (mZ-g[iPop]+1):mZ, 3, replace=FALSE )
			z <- sapply( 1:2, function(zi){Z[,rrGen[zi],rrcPar[iGenThin,iChain,zi] ]})
			tmp.i = tmp.i + 1
		}
		dz <- (z[,1]-z[,2])	#jump vector as the difference between two random states (from z2 towards z1)
		if( !is.null(ctrl$probUpDir) && (ctrl$probUpDir != 1/2) ){
			lz <- sapply( 1:2, function(zi){ rLogDen[ rrGen[zi], rrcPar[iGenThin,iChain,zi] ] })#get the logDensity for the random states
			if( all(is.finite(lz)) && (lz[1] < lz[2]) )	# when proposing a jump between two states towards lower density
				if( runif(1) >= (2*ctrl$probUpDir-1) ) dz = -dz	# increase chance of going directio of upward LogDensity
		}
		if ( runif(1)< ctrl$pGamma1 ) { 
			gamma_par = ctrl$F1 # to be able to jump between modes
		} else {
			gamma_par = ctrl$F2 * runif(d$parms, min=1-ctrl$epsMult, max=1+ctrl$epsMult)    # multiplicative error to be applied to the difference
			# gamma_par = F2 
		}
		if (ctrl$epsAdd ==0) {  # avoid generating normal random variates if possible
			#xPropChain = X[,iChain] + gamma_par * dz 
			xStepChain = gamma_par * dz 
		} else {
			#xPropChain = X[,iChain] + gamma_par * dz  + rnorm(d$parms,0,ctrl$epsAdd)
			xStepChain = gamma_par * dz  + rnorm(d$parms,0,ctrl$epsAdd)
		}
		rExtra=0
	} #DE-parallel direction update
	#c( xPropChain, rExtra )
	c( xStepChain, rExtra )
	### vector with first components the proposal and last component some extra LogDensity for decision for snooker decision
}


.doDEMCStep <- function( 
	### version before sfRemoteWrapper, Perfrom one DEMC step, function to be  alled in remote process
	x, logDenCompX, logDenX, 
	step, rExtra, 
	temp,			##<< temperature of the current step and population
	argsDEMCStep	
### list with components \describe{
### \item{fDiscrProp,argsFDiscrProp}{function and additional arguments applied to xProp, e.g. to round it to discrete values}
### \item{argsFLogDen, fLogDenScale}{additional arguments to fLogDen and scalar factor applied to result of fLogDen}
### \item{}{}
### }
){
	##details<< 
	## If argsDEMCStep is.name, it is evaluated in remote process. This allows to distribute parameters only once
	## per \code{sfExport("argsDEMCStep")} and calling this function with \code{argsDEMCStep=as.name("argsDEMCStep")}
	if( is.name(argsDEMCStep)) argsDEMCStep=eval.parent(argsDEMCStep)	#do not need to be passed each time
	boResFLogDenX <- { .nc <- ncol(logDenCompX); ( !is.null(.nc) && (.nc > 0) ) }	#if number of columns > 0
	body <- expression( with( argsDEMCStep, {
				accepted<-FALSE
				xProp = x + step
				if( is.function(fDiscrProp)) # possiblity to discretize proposal
					xProp = do.call(fDiscrProp,xProp,argsFDiscrProp, quote=TRUE)
				res <- if(boResFLogDenX)
						do.call( fLogDen, c(list(xProp), list(logDenCompX), argsFLogDen) )	# evaluate logDen
					else
						do.call( fLogDen, c(list(xProp), argsFLogDen) )	# evaluate logDen
				logDenProp <- sum(res)*fLogDenScale		# maybe a vector, reduce to scalar
				if( is.finite(logDenProp) ){
					logDenCompProp <- if( 0 < length(logDenCompX) )	#extract the components of res that are handled internally
						if( !all( names(logDenCompX) %in% names(res) )){ 
							stop(paste("not all names of logDenCompX in return of fLogDen: ",paste(names(res),collapse=",")))
							res[names(logDenCompX)]
						}else FLogDenX
					logr = (logDenProp+rExtra - logDenX) / temp
					#Metropolis step
					if ( is.finite(logr) & (logr) > log(runif(1)) ){
						accepted <- TRUE
						x <- xProp
						logDenCompX <- logDenCompProp
						logDenX <- logDenProp
					}
				}
				list(accepted=accepted,x=x,logDenCompX=logDenCompX,logDenX=logDenX)
			})) 
	##details<<  
	## if argsDEMCStep entry remoteDumpfileBasename is characters
	## a dump is created in this file when an error occurs and execution stops with an error
	## see ?dump.frames
	res <- if( !is.null(argsDEMCStep$remoteDumpfileBasename)){
			res <- try( eval(body) )
			if( inherits(res, "try-error") ){
				dump.frames(argsDEMCStep$remoteDumpfileBasename,TRUE)
				stop(res)
			}else res
		}else{
			eval(body)
		}
	### list with components \describe{
	### \item{accepted}{boolean scalar: if step was accepted}
	### \item{x}{numeric vector: current position in parameter space}
	### \item{logDenCompX}{numeric vector: result of fLogDen for current position}
	### \item{logDenX}{numeric scalar: logDen of current position}}
}

