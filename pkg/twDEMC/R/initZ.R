initZtwDEMCNormal <- function(
	### Generate an initial population of states for \code{\link{twDEMCInt}}.
	thetaPrior,	
	### vector of parameters, point estimate
	covarTheta, 
	### the a prior covariance of parameters 
	nChains=8, 
	### number of chains to run
	nPops=2, 
	### number of independent populations among the chains 
	m0=calcM0twDEMC(length(thetaPrior),nPops,nChains)
	### number of initial states for each chain
	,doIncludePrior=TRUE
	### If TRUE, then last sample of chain 1 to the prior estimate, 
	### which might be already a kind of best estimates by an optimization. 
){
	# InittwDEMC

	##seealso<<  
	## \code{\link{twDEMCInt}}
	## \code{\link{calcM0twDEMC}}
	
	##details<< 
	## There are several methods to establish the initial population for a twDEMC run. \itemize{
	## \item{ drawing from a multivariate normal distribution: this method  } 
	## \item{ subsetting the result of a former twDEMC run: \code{\link{initZtwDEMCSub.twDEMC}}  } 
	## \item{ extending the result of a former twDEMC run to include more parameters: \code{\link{initZtwDEMCExt.twDEMC}}  } 
	## \item{ selecting the N closes points from a sequence of points in parameter space \code{\link{constrainNStack}}  } 
	## \item{ selecting the points inside a confindenc ellipsis in parameter \code{\link{constrainCfStack}}  } 
	##}
	
	#?rmvnorm #in package mvtnorm
	if( 0==length(thetaPrior))
		stop("initZtwDEMCNormal: no parameters given")
	if( 0==length(names(thetaPrior))){ 
		warning("initZtwDEMCNormal: no names of parameters provided")
		names(thetaPrior) <- paste("theta",1:length(thetaPrior),sep="") 
	}
	Zinit <- if( length(thetaPrior)==1 ){
		abind( lapply( 1:nChains, function(i){ t(matrix(rnorm( m0, mean=thetaPrior, sd=covarTheta), dimnames=list(NULL,names(thetaPrior)) )) }), along=3 )
	}else{
		abind( lapply( 1:nChains, function(i){ t(rmvnorm( m0, mean=thetaPrior, sigma=covarTheta )) }), along=3 )
	}
	# Set the last sample of chain 1 to the prior estimate, which might be already a kind of best estimates by an optimization.
	if( doIncludePrior )
		Zinit[,m0,1 ] <- thetaPrior
	dimnames(Zinit) = list(parms=names(thetaPrior),steps=NULL,chains=NULL)
	Zinit
	### a matrix of number of parameters by number of individuals (d x m0 x Npop), with d dimension of theta
}

calcM0twDEMC <- function(
	### Calculate appropriate number of cases for initializing twDEMC.
	nPar,	##<< the number of parameters to estimate
	nPops, 	##<< the number of independent populations
	nChains ##<< the number of chains 
){
	##seealso<<   
	## \code{\link{initZtwDEMCNormal}}
	
	##details<< see terBraak 2006 and 2008
	res <- max(4,(8*nPar*nPops)%/%nChains)
	res
	### length of each chain so that each population is initialized with 8*nPar cases 
}

setMethodS3("initZtwDEMCSub","matrix", function(
	### generates an appropriate initial sample of parameter vectors for twDEMC from subsampling an array
	Zinit1,					##<< the mcmc matrix to subsample (column variable, rows cases) 
	vars=colnames(Zinit1),	##<< which variables to keep
	nChains=4,	
	nPops=1, 
	m0=calcM0twDEMC(length(unique(vars)),nPops,nChains), ##<< number of required cases for initialization
	... 
){
	##seealso<<   
	## \code{\link{initZtwDEMCNormal}}
	rrc <- lapply(1:nChains, function(iChain){ sample( (1:nrow(Zinit1)), m0, replace=TRUE)} )
	res <- abind( lapply( rrc, function(rr){ t(Zinit1[rr,vars,drop=FALSE])} ),rev.along=0 )
	dimnames(res) <- list( parms=vars, steps=NULL, chains=NULL ) 
	res
	### an array of dimension suitable for Zinit for twDEMCInt
})

setMethodS3("initZtwDEMCSub","twDEMC", function(
		### generates an appropriate initial sample of parameter vectors for twDEMC from subsampling a previous result
	vtwdemc,	##<< the twDEMC list to subsample 
	vars=rownames(vtwdemc$parms),	##<< which variables to keep
	nChains=ncol(vtwdemc$rLogLik),	
	nPops=ncol(vtwdemc$temp),
	...	## other parameters passed to initZtwDEMCSub.default, e.g. m0 
){
	##seealso<<   
	## \code{\link{initZtwDEMCNormal}}
	Zinit1 <- stackChains(vtwdemc)[,-1]	#one big chain
	initZtwDEMCSub.matrix(Zinit1,vars,nChains,nPops,...)
})


setMethodS3("initZtwDEMCExt","matrix", function( 
	### subsampling and extending twDEMC with new variables
	Zinit1,	##<< the twDEMC list to subsample 
	thetaPrior,	##<< numeric vector: mu of multivariate gaussian prior distribtuion 
	covarTheta,		##<< numeric vector: sigma of multivariate gaussian prior distrbituion 
	nChains=4,	
	nPops=1, 
	m0=calcM0twDEMC(length(thetaPrior),nPops,nChains),
		### number of required cases 
	...
){
	# initZtwDEMCExt
	##details<< 
	## the variables in thetaPrior that are part of vtwdemc are subsampled
	## the other variables are drawn from prior distribution
	## assuming no correlations between variables present and absent in vtwdemc and 
	if( is.null(names(thetaPrior)) )
		stop("initZtwDEMCExtDEMCzsp: thetaPrior provided without names attribute")
	#mapping index thetaPrior to column index in pss1
	pssInd <- as.numeric(lapply( 1:length(thetaPrior), function(i){ which( names(thetaPrior)[i] == colnames(Zinit1))}))
	mcPars <- which( !is.na(pssInd) )   
	extPars <- (1:length(thetaPrior))[ -mcPars ]
	Zinit <- array(NA, dim=c(length(thetaPrior), m0, nChains) )
	if( length(mcPars) > 0 ){
		ZinitMc <- abind( lapply( 1:nChains, function(i){ t( Zinit1[sample(1:nrow(Zinit1), m0, replace=TRUE), pssInd[mcPars], drop=FALSE ]) }), along=3 )
		Zinit[mcPars,,] <-  ZinitMc
	}		
	if( length(extPars) > 0){
		ZinitZtwDEMCExt <- initZtwDEMCNormal( thetaPrior[extPars], covarTheta[extPars,extPars,drop=FALSE], nChains=nChains, nPops=nPops,m0=m0 )
		Zinit[extPars,,] <- ZinitZtwDEMCExt
	}
	dimnames(Zinit) <- list( parms=names(thetaPrior), steps=NULL, chains=NULL ) 
	Zinit
	
	##seealso<<   
	## \code{\link{initZtwDEMCNormal}}
})

setMethodS3("initZtwDEMCExt","twDEMC", function( 
	 ### subsampling and extending twDEMC with new variables
	vtwdemc,	##<< the twDEMC list to subsample 
	thetaPrior,	##<< numeric vector: mu of multivariate gaussian prior distribtuion 
	covarTheta,		##<< numeric vector: sigma of multivariate gaussian prior distrbituion 
	nChains=ncol(vtwdemc$rLogLik),	
	nPops=ncol(vtwdemc$temp),
	...		##<< e.g. m0
){
	##seealso<<   
	## \code{\link{initZtwDEMCNormal}}
	Zinit1 <- stackChains(vtwdemc)[,-1]	#one big chain
	initZtwDEMCExt.matrix( Zinit1, thetaPrior, covarTheta, nChains=nChains, nPops=nPops, ... )
})


constrainNStack <- function( 
	### Subsetting a sequence of parameter vectors. Keeps only the n cases in pss1 that are closest to thetaPrior.
	pss1, 					##<< numeric matrix: the stack to constrain, see details 
	thetaPrior,				##<< the point in parameter spcae for which to select the closest values 
	n=nrow(pss1)%/%4, 		##<< the number of rows in output, defaults to 1/4 of nrow(pss1)
	vars=names(thetaPrior), ##<< names or indices to constrain thetaPrior and invCovarTheta
	invCovarTheta=if( 0<length(thetaPrior)) solve(cov(pss1[,names(thetaPrior[vars]),drop=FALSE])) else numeric(0),
		### the inverse of the covaraince matrix for thetaPrior, defaults to inverse of sample covariance
	returnAlpha=FALSE 		##<< switch to return also significance level
){
	# constrainNStack

	##seealso<<   
	## \code{\link{initZtwDEMCNormal}}
	## \code{\link{constrainCfStack}}
	
	##details<< 
	## pss1 is a matrix with columns corresponding to variables 
	## and rows corresponding to cases. 
	## It is typicalla the result of \code{\link{stackChains.twDEMC}}
	if( !(n < nrow(pss1)) ){
		warning("constrainNStack: number of rows in pss1 is not greater than n, returning full pss1")
		return( pss1)
	}
	tmp.ind <- if( 0<length(thetaPrior)){
		if( !(all(names(thetaPrior[vars]) %in% colnames(pss1))) ){
			stop("constrainNStack: all components of thetaPrior[vars] must be in colnames(pss1)")
		}
		.invCovarTheta <- if( is.null(colnames(invCovarTheta)) | is.null(rownames(invCovarTheta)) ){
			.lv <- length(vars)
			if( !identical( c(.lv,.lv), dim(invCovarTheta)) )
				stop("supplied invCovarTheta without dimensions and nonmatching dimensions")
			invCovarTheta
		}else{
			if( !all(vars %in% colnames(invCovarTheta)) | !all(vars %in% rownames(invCovarTheta)) )
				stop("supplied invCovarTheta with dimnames not including all names of thetaPrior")
			invCovarTheta[vars,vars]
		}
		# calculate the distance and order by that
		tmp.d <- sapply(vars, function(var){ pss1[,var] - (thetaPrior[var])})
		#plot( tmp.d[,1], tmp.d[,2] )
		if( length(vars) > 1) 
			tmp.mahalanobis <- apply( tmp.d, 1, function(tmp.d){ t(tmp.d) %*% .invCovarTheta %*% tmp.d })
		else 
			tmp.mahalanobis = tmp.d^2 * 1/var(pss1[,vars]) 
		tmp.ind <- order(tmp.mahalanobis)[1:n]
	}else{
		# no thetaPrior given: just subsample
		tmp.ind <- sample.int(nrow(pss1),n) 
	}
	psc <- pss1[ tmp.ind, ]
	if( returnAlpha){
		tmp.alpha <- if( 0 < length(thetaPrior)){
			# calculate the corresponding alpha (p) of the confidence ellipsis
			# see http://www.stat.psu.edu/online/development/stat505/05_multnorm/04_multnorm_geom.html
			tmp.maxdist = tmp.mahalanobis[tmp.ind[n] ]
			tmp.alpha = pchisq( q=tmp.maxdist, df=length(vars) )
		}else 1
		list( res=psc, alpha=tmp.alpha)
	}else{
		psc
	}
	### pss1[closestValues,]
	### if returnAlpha=TRUE then a list is returned with \describe{ 
	### \item{res}{ value as above }
	### \item{alpha}{ the significance level of the corresponding confidence ellipsis }} 
}
#mtrace(constrainNStack)
#twUtestF(constrainNStack,test="test.allVars")

constrainCfStack <- function( 
	### Subsetting a sequence of parameter vectors. Keeps only the the cases in pss1 that are inside a confidence ellipsis around thetaPrior.
	pss1, 			##<< numeric matrix: the stack to constrain, see details 
	thetaPrior,	##<< the point in parameter spcae for which to select the closest values 
	alpha=0.95, ##<< the conficence range of the ellipsis 
	vars=names(thetaPrior), ##<< names or indices to constrain thetaPrior and invCovarTheta
	invCovarTheta=if(0<length(thetaPrior)) solve(cov(pss1[,vars,drop=FALSE])) else numeric(0)
	### the inverse of the covaraince matrix for thetaPrior, defaults to inverse of sample covariance
){
	# constrainCfStack
	
	##seealso<<   
	## \code{\link{initZtwDEMCNormal}}
	## \code{\link{constrainNStack}}
	if( 0==length(thetaPrior)) return( pss1 )	# if no thetaPrior is given, just return the original matrix 
	if( !(all(names(thetaPrior[vars]) %in% colnames(pss1))) ){
		stop("constrainNStack: all components of thetaPrior[vars] must be in colnames(pss1)")
	}
	.invCovarTheta <- if( is.null(colnames(invCovarTheta)) | is.null(rownames(invCovarTheta)) ){
			.lv <- length(vars)
			if( !identical( c(.lv,.lv), dim(invCovarTheta)) )
				stop("supplied invCovarTheta without dimensions and nonmatching dimensions")
			invCovarTheta
		}else{
			if( !all(vars %in% colnames(invCovarTheta)) | !all(vars %in% rownames(invCovarTheta)) )
				stop("supplied invCovarTheta with dimnames not including all names of thetaPrior")
			invCovarTheta[vars,vars,drop=FALSE]
		}
	#first remove all values outside the 95% confidence ellipsis
	# see http://www.stat.psu.edu/online/development/stat505/05_multnorm/04_multnorm_geom.html
	tmp.d <- sapply(vars, function(var){ pss1[,var] - (thetaPrior[var])})
	tmp.mahalanobis <- apply( tmp.d, 1, function(tmp.d){ t(tmp.d) %*% .invCovarTheta %*% tmp.d })
	tmp.chisq = qchisq( p=alpha, df=length(vars) )
	psc <- pss1[ tmp.mahalanobis <= tmp.chisq, ]
	if( nrow(psc) < 10 )
		stop(paste("constrained too strong, selected nRows=",nrow(psc)))
	psc
	### pss1[ withinConfidenceInterval, ]	
}


replaceZinitCases <- function( 
	### Replaces states of Zinit that yield non-finite rLogLik by sampling other states.
	Zinit, ##<< initial states of the form required by \code{\link{twDEMCInt}} 
	boMat ##<< boolean matrix with rows cases and columns parameters
){
	##seealso<<   
	## \code{\link{initZtwDEMCNormal}}
	## \code{\link{replaceZinitNonFiniteLogLiks}}
	
	##details<< 
	## Samples for the first half of chains are sampled from good cases of the second half of chains.
	## Samples for the second half of chains are sampled from good cases of first half of chains.
	tmp.nChains <- dim(Zinit)[3]
	if( !identical(dim(Zinit)[c(2,3)], dim(boMat)) )
		stop("dimensions of Zinti and boMat do not match")
	if( !is.logical(boMat) )
		stop("boMat is not a logical matrix")
	goodCases <- which( boMat, arr.ind=TRUE )
	goodCases1 <- goodCases[ goodCases[,2] <= tmp.nChains/2, ,drop=FALSE]	#non-missing values of first half
	goodCases2 <- goodCases[ goodCases[,2] > tmp.nChains/2, ,drop=FALSE]	#non-missing values of second half
	badCases <- which( !boMat, arr.ind=TRUE )
	badCases1 <- badCases[ badCases[,2] <= tmp.nChains/2, ,drop=FALSE] #missing values of the first half
	badCases2 <- badCases[ badCases[,2] > tmp.nChains/2, ,drop=FALSE] #missing values of the second half
	#replace non-finite yielding states of first halv by sampling the good cases rows of the second half
	tmp.j = sample( 1:nrow(goodCases2), nrow(badCases1), replace=TRUE )
	for( i in seq(along.with=badCases1[,1]) )	Zinit[, badCases1[i,1], badCases1[i,2] ] = Zinit[, goodCases2[tmp.j[i],1], goodCases2[tmp.j[i],2]  ]
	#replace non-finite yielding states of second halv by sampling the good cases rows of the first half
	tmp.j1 = sample( 1:nrow(goodCases1), nrow(badCases2), replace=TRUE )
	for( i in seq(along.with=badCases2[,1]) )	Zinit[, badCases2[i,1], badCases2[i,2] ] = Zinit[, goodCases1[tmp.j1[i],1], goodCases1[tmp.j1[i],2]  ]
	Zinit
	### Zinit, with several cols (parameter vectors) replaced by other cols
}

replaceZinitNonFiniteLogLiks <- function( 
	### Replaces states of Zinit that yield non-finite rLogLik by sampling other states.
	Zinit, ##<< initial states see InitDEMCzsp 
	rLogLik ##<< tmp.rLogLik: calculated logLikelihoods for all the states in Zinit (rows cases and columns parameters) 
){
	##seealso<<   
	## \code{\link{initZtwDEMCNormal}}
	## \code{\link{replaceZinitCases}}
	
	##details<< 
	## In order for twDEMC to start up effectively, it is important that chains start from values, where the logLikelihood is finite
	bo <- is.finite(rLogLik)
	dim(bo) <- dim(Zinit)[2:3]
	replaceZinitCases( Zinit, bo )
	### Zinit, with several cols (parameter vectors) replaced by other cols
}
#twUtestF(replaceZinitNonFiniteLogLiks)
attr(replaceZinitNonFiniteLogLiks,"ex") <- function(){
	data(twLinreg1)
	attach( twLinreg1 )
	mtrace(initZtwDEMCNormal)
	Zinit <- initZtwDEMCNormal( theta0, diag(4*sdTheta^2), nChains=8, nPops=2)
	dim(Zinit)
	res <- twCalcLogLikPar(logLikGaussian, stackChains(Zinit)
		,fModel=dummyTwDEMCModel		### the model function, which predicts the output based on theta 
		,obs=obs			### vector of data to compare with
		,invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)$logLik
	plot(density(res))
	res[res< -30] <- NA
	resM <- matrix(res, ncol=dim(twdemcEx1$parms)[3])
	Zinit2 <- replaceZinitNonFiniteLogLiks(Zinit, resM)
}

replaceZinitNonFiniteLogLiksLastStep <- function( 
	### Replaces states of last step, i.e. column of Zinit that yield non-finite rLogLik by sampling other states.
	Zinit 			##<< initial states see InitDEMCzsp
	,fLogLik 		##<< the logLik Function
	,nPops=1		##<< number of populations. States are only choosen from same population
	,iStep=dim(Zinit)[2]	##<< the step for which to replace nonfinite yielding parameters.
	,maxSteps=16
	,...			##<< arguments to \code{\link{twCalcLogLikPar}}
){
	##seealso<<   
	## \code{\link{initZtwDEMCNormal}}
	res <- twCalcLogLikPar( fLogLik, t(adrop(Zinit[,iStep,,drop=FALSE],2)), ... )
	rLogLik <- res$logLik
	rLogLikComp <- res$logLikComp
	resFLogLik <- res$resFLogLik
	iNonfinite <- which( !is.finite(rLogLik) )
	nReplace <- length(iNonfinite)
	if( nReplace>0 ){ 
		# matrix finLokLik markes the states which we do not want to sample again by FALSE
		dimZinit <- dim(Zinit)
		finLogLik <- matrix( TRUE, nrow=dimZinit[2], ncol=dimZinit[3] )
		finLogLik[iStep,] <- FALSE	#do not choose replacements from this step
		chainsPops <- matrix( 1:dimZinit[3], ncol=nPops )
		nChainsPop <- nrow(chainsPops)
		rChain <- integer(nReplace)
		rStep <- integer(nReplace)
		i=1
		while( (nReplace>0) & i<=maxSteps ){
			for( j in 1:nReplace){
				chainsPop <- chainsPops[,((iNonfinite[j]-1)%/%nChainsPop)+1]
				rChain[j] <- sample(  chainsPop, 1 )
				stepsPop <- which( finLogLik[,rChain[j] ] )
				rStep[j] <- sample(stepsPop , 1 )
				Zinit[,iStep, iNonfinite[j] ] <- Zinit[,rStep[j],rChain[j] ]
				finLogLik[rStep[j],rChain[j] ] <- FALSE		#do not sample those states again
			}
			res <- twCalcLogLikPar( fLogLik, t(adrop(Zinit[,iStep,iNonfinite,drop=FALSE],2)), ... ) 
			rLogLik[iNonfinite] <- res$logLik
			rLogLikComp[iNonfinite,] <- res$logLikComp
			if( 0 < length(resFLogLik))
				resFLogLik[iNonfinite,] <- res$resFLogLik
			iNonfinite <- which( !is.finite(rLogLik) )
			nReplace <- length(iNonfinite)
			i=i+1
		}
		if( nReplace > 0 )
			warning("could not replace all non-finite states in row.")
	}
	list(Zinit=Zinit, rLogLik=rLogLik, rLogLikComp=rLogLikComp, resFLogLik=resFLogLik)
	### list with components \describe{
	### \item{Zinit}{given Zinit with some states in last row replaced by other states from Zinit.}
	### \item{rLogLik, rLogLikComp, resFLogLik}{numeric matrix: results of \code{\link{twCalcLogLikPar}} for last row for all chains}
	### }
}
#twUtestF(replaceZinitNonFiniteLogLiks)







