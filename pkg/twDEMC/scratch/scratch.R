#(tres <- twUtestF(combinePops,"test.twDEMC.apply"))
twDEMCApply <- function(
	### applies function to each chain in $parms of an twDEMC
	vtwdemc,
	fun,
	### function(x,name) applied to each variable in twDEMC
	otherMargin = 1,
	### margin for that each index the an array slice is provided as x to fun, defaults to 1: Vars (2: generations, 3:chains)
	...
){
	##details<< 
	## the arguments of fun are \describe{
	## \item{x}{a matrix corresponging to margin (default with rows correponding to generations and cols corresponding to) }
	## \item{name}{the name of the remaining margin, or the number of the item for otherMargin=2,3}}
	## the return value must be of same size as x
	if( length(otherMargin)!= 1) 
		stop(paste("twDEMC.applyParms: margin must be of type c(i,j) with i!=j and i,j in 1:3"))
	res <- vtwdemc
	#for( iVarName in unlist(dimnames(vtwdemc$parms)[otherMargin]) ){
	.ind <- 1:(dim(vtwdemc$parms)[otherMargin])
	.names <- unlist(dimnames(vtwdemc$parms)[otherMargin])
	if( is.null(.names) ) .names=.ind
	for( iVar in .ind ){
		#aperm(res$parms,c(otherMargin,margin)) <- fun( aperm(vtwdemc$parms,c(otherMargin,margin))[iVar,,],iVar,...)
		iVarName = .names[iVar]		
		if( otherMargin==1) res$parms[iVar,,] <- fun( vtwdemc$parms[iVar,,],iVarName,...)
		else if( otherMargin==2) res$parms[,iVar,] <- fun( vtwdemc$parms[,iVar,],iVarName,...)
		else if( otherMargin==3) res$parms[,,iVar] <- fun( vtwdemc$parms[,,iVar],iVarName,...)
	}
	res
}
#mtrace(twDEMCApply)

test.twDEMCApply <- function(){
	.varDistr <- c(a="lognorm",b=NULL)
	.tmp <- twDEMCApply(twdemcEx1,function(x,varname){ transOrigPopt.default(x,.varDistr[varname])})
	checkEquals(dim(twdemcEx1$parms),dim(.tmp$parms))
	checkEquals(.tmp$parms["a",,], exp(twdemcEx1$parms["a",,]) )
	checkEquals(.tmp$parms["b",,], (twdemcEx1$parms["b",,]) )
	
	.tmp2 <- twDEMCApply(twdemcEx1,otherMargin=2, fun=function(x,varname){ x+1 })
	checkEquals(.tmp2$parms[,,], (twdemcEx1$parms[,,]+1) )
}

twDEMC <- function(
	### Differential Evolution Markov Chain 
	Zinit,...
){	
	##details<< 
	## calls correct initializer of \code{\link{twDEMCBlockInt}} depending on the class of Zinit
	UseMethod("twDEMC"
	##seealso<< 
	## \code{\link{twDEMCBlockInt}}
	## \code{\link{twDEMC.default}}
	## \code{\link{twDEMC.twDEMC}}
	)}

#gamma <- popt.tmp; gamma0 <- tmp.m
rotateInd <- function(gamma, gamma0, A){
	# tansInd
	#
	# transpose parameter vector gamma to independent parameters xi
	# gamma parameters to transpose
	# gamma0 approximate mean of the posterior
	# A matrix square root, usually obtained by chol(cov), where cov is the covariance of the posterior
	gd <- gamma - gamma0
	xi <- backsolve( A, gd )
	xi
}
rotateIndBack <- function(xi,gamma0,A){
	# tansIndBack
	#
	# transpose independent parameters xi back to original parameters (inverse of rotateInd)
	# xi parameters to transpose back
	# gamma0 approximate mean of the posterior
	# A upper diagonal matrix square root, usually obtained by chol(cov), where cov is the covariance of the posterior
	gamma = A %*% xi + gamma0
	gamma
}

#------------ rotating parameter vectors
tmp.f <- function(){
	?chol	
	tmp.cov <- cov( pss[,-1] )	#sample covariance
	tmp.m <- colMeans(pss[,-1])	#sample mean
	A = chol(tmp.cov)		#upper t
	
	xi <- transInd(popt.tmp,tmp.m, A)
	xi2 <- solve(A) %*% gd
	xi - xi2
	tmp <- rotateIndBack(xi,tmp.m, A)
	popt.tmp - tmp
	
	mtrace(transInd)
	pssXi <- t(transInd( t(pss[sample(1:nrow(pss),20),-1]), tmp.m, A ))
	cor(pssXi)
	corOrdered(pssXi)	#why still highly dependent?
}



checkConvergenceFalse <- function(res, addArgs=list() ){ FALSE }

checkConvergenceGelman <- function(res, addArgs=list() ){
	# checkConvergenceGelman
	#
	# gelman RHat criterion applied to second half of each chain of a population
	# res see return value of DEMCzsp ($parms (d x nStep x nChain) )
	# value: TRUE/FALSE
	# see Gelman04 (twutz:Gelman04_3#Inference_and_assessing_convergence)
	# addArgs
	# $burninFrac fraction of the chain to be discarded (default 0.5)
	# $rHatMin rHat criterion, returns TRUE if smaller than this 
	if( is.null(addArgs$burninFrac) ) addArgs$burninFrac=0.5
	if( is.null(addArgs$rHatMin) ) addArgs$rHatMin = 1.1
	l <- dim(res$parms)[2]
	res2 <- res$parms[ ,(l*addArgs$burninFrac):l, ]	# second half of all the chains
	
	n <- dim(res2)[2]	# number of steps
	m <- dim(res2)[3]   # number of chains
	rl <- sapply( 1:dim(res2)[1], function(vn){	#over all variables (estimands)
			muyj <- sapply(1:m, function(j){ mean(res2[vn,,j]) }) # mean of each chain
			muy <- mean(muyj)	# mean across chains
			B <- n/(m-1)*sum((muyj-muy)^2) # between chain variance
			s2j <- sapply(1:m, function(j){sum((res2[vn,,j]-muyj[j])^2) })/(n-1) #squared sums of each chain
			W <- sum(s2j)/m	# within chain variance
			VarPlus=(n-1)/n*W + 1/n*B
			#c( B=B, W=W, VarPlus=VarPlus, Rhat=round(sqrt(VarPlus/W),3), effSize=min(n*m,n*m*VarPlus/B) )
			c(  Rhat2=VarPlus/W )
		})
	all(rl <= (addArgs$rHatMin)^2 )
}
#mtrace( checkConvergenceGelman )

.testVectorApply <- function(){
	tmp.f <- function(x){ 1:4 }
	x <- matrix(1:3, nrow=3, ncol=5)
	(.tmp <- apply( x, 2, tmp.f ))	#each vector becomes a column
	
	tmp.f1 <- function(x){ 1 }
	(.tmp <- apply( x, 2, tmp.f1 ))	#each scaler becomes component of a vector
	
}

.tmp.f <- function(
	doLocalAdaptation=FALSE
	### if TRUE, then states for differential proposal are sampled with heigher weight from current chain and recent past
	, doAdaptT=FALSE
### if TRUE, Temperature will be adapted to acceptance rate
){	
}


.doDEMCSteps <- function( iGenT0, ctrl, 
	X, xStep, 
	fDiscrProp, argsFDiscrProp, 
	logDenCompX, argsFLogDen, fLogDenScale, debugSequential, ...
){
	for( iGenT in (1:ctrl$thin) ){
		iGen = iGenT0+iGenT
		xProp = X + xStep[,,iGenT]
		# if necessary, discretize proposal
		if( is.function(fDiscrProp))
			xProp = do.call(fDiscrProp,xProp,argsFDiscrProp, quote=TRUE)
		
		##details<< 
		## in order to support a two-level Metropolis desition
		#twCalcLogDenPar expects parameters in columns, need transpose
		.resLogDenPar <- twCalcLogDenPar(fLogDen=fLogDen, xProp=t(xProp), logDenCompX=logDenCompX, argsFLogDen=argsFLogDen, fLogDenScale=fLogDenScale, debugSequential=debugSequential, ...)
		logfitness_x_prop <- .resLogDenPar$logDen
		logDenCompProp <- .resLogDenPar$logDenComp
		
		TcurStep = TstepFixed[iGen,]
		#logr =  (logfitness_x_prop+rExtra - logfitness_X)/Tstep[iGen]
		logr =  (logfitness_x_prop+rExtra[,iGenT] - logfitness_X) / rep(TcurStep, each=nChainPop)
		# print(c(logfitness_X[i], logfitness_x_prop ,logr,rExtra))
		
		#Metropolis step for each chain
		for (i in iseq){	#body must affect outer variables: cannot use apply instead of for
			if ( is.finite(logr[i]) & (logr[i]) > log(runif(1)) ){
				acceptN[i] = acceptN[i]+1
				X[,i] = xProp[,i]
				logfitness_X[i] = logfitness_x_prop[i]
				if( 0<length(logDenCompX) ) 
					logDenCompX[i,] = logDenCompProp[i,]	#here chains are rows
				acceptWindow[ acceptPos, i ] <- TRUE
			}
		} #Metropolis step for each chain
	}# iGenT within thinning interval
	resDo <- list(	acceptN=acceptN, X, logfitness_X, logDenCompX, acceptWindow )
	### list with components \describe{
	### \item{acceptN}{vector number of accepted steps per chain in thinning interval}
	### \item{X}{matrix current position for each chain (column?)}
	### \item{logfitness_X}{vector current logDen of chains}
	### \item{logDenCompX}{matrix: result of fLogDen for last accepted state per chain}
	### \item{acceptWindow}{}
}

.generateXPropThin <- function(
	nPops,
	Z,mZ,
	ctrl,
	...
){
	d <- as.list(structure(dim(Z),names=c("parms","gen","chains")))
	d$gen <- mZ
	d$steps <- ctrl$thin
	nChainPop = d$chains %/% nPops
	##details<<  
	## Random states for chains for difference vector generation are within subsets of populations.
	## This allows simulating independent population of chains.
	## The acceptance rate may differ amonng populations. Hence, the set of previous generations to 
	## randomly select from also differs between poplations.
	# integer array (thinSteps*nChain*3) sample of chains, within populations
	rrChains <- abind( lapply( 1:nPops, function(iPop){
				sChains <- ((iPop-1)*nChainPop+1):(iPop*nChainPop)
				rrChainsPop <- array( sample(sChains, d$steps*nChainPop*3, replace = TRUE), dim=list(gen=d$steps,chain=nChainPop,iPop=3) )
			}), along=2 )
	# integer array same dimension of sample of generations
	rrGen <- abind( lapply( 1:nPops, function(iPop){
				sGens <- (mZ-g[iPop]+1):mZ 
				rrGenPop <- array( sample(sGens, d$steps*nChainPop*3, replace = TRUE), dim=list(gen=d$steps,chain=nChainPop,i=3) )
			}), along=2 )
}

calcCorrParsDf <- function(parms, normpopt){
	# calcCorrParsDf
	#
	# caluclates the values of correlated parameters on transformed, i.e. normal, scale
	# based on dataframe parms$corrParsDf, which has columns y,x (char), a, b
	# normpcorr[y] = a*normpopt[x] +b
	# which is the result of linApprox in R/sensitivity.R of this package
	tmp5 <- parms$corrParsDf
	if( is.data.frame(tmp5) && (nrow(tmp5) > 0)){
		tmp <- sapply( 1:nrow(tmp5), function(i){
				row = tmp5[i,]
				row$a * normpopt[ row$x ] + row$b 
			})
		names(tmp) <- tmp5$y
		tmp
	}else c()
}

calcCorrParsDfMat <- function(parms, normpopt){
	# calcCorrParsDfMat
	#
	# same as calcCorrParsDf but with normpopt specified as a matrix with mutliple rows
	# caluclates the values of correlated parameters on transformed, i.e. normal, scale
	# based on dataframe parms$corrParsDf, which has columns y,x (char), a, b
	# normpcorr[y] = a*normpopt[x] +b
	# which is the result of linApprox in R/sensitivity.R of this package
	tmp5 <- parms$corrParsDf
	if( is.data.frame(tmp5) && (nrow(tmp5) > 0)){
		tmp <- sapply( 1:nrow(tmp5), function(i){
				row = tmp5[i,]
				row$a * normpopt[, row$x ] + row$b 
			})
		colnames(tmp) <- tmp5$y
		tmp
	}else c()
}



.test.ggstat <- function(){
	str(diamonds)
	p1 <- ggplot(diamonds, aes(x=color,y=carat))
	p2 <- p1 + geom_point()
	p1 + stat_summary(aes(colour="mean"), fun.y=mean, geom="point")
	+ scale_colour_hue("Mean")
	
	x <- rnorm(100) 
	mtrace(dnorm)
	p1 <- ggplot(x)
	p1 + geom_
	qplot(x, geom="density") + stat_function(fun = dnorm, colour="red") 	
}

setMethodS3("transOrigPopt","twDEMC", function( 
		### Applies \code{\link{transOrigPopt.default}} to each column of parameters in \code{vtwdemc}.
		vtwdemc, 
		### list of class twDEMC with $parms in transformed scale 
		poptDistr=eval(parse(text="parDistr$trans")),
		### character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt
		...
	){
		##seealso<<   
		## \code{\link{transOrigPopt.default}}
		res <- vtwdemc
		normpopt <- vtwdemc$parms
		##details<< 
		## Either \code{poptDistr} has names for each column name,
		## or \code{poptDistr} has the same length as \code{colnames(normpopt)}.
		if( !is.null(names(poptDistr)) )
			if( all(rownames(normpopt) %in% names(poptDistr)) )
				poptDistr=poptDistr[ rownames(normpopt) ]
		if( length(poptDistr) != nrow(normpopt) )
			stop("poptDistr must have entries for each row or must have the same length as nrow(normpopt$parms).")
		for( vi in seq(along.with=rownames(normpopt)) ){
			res$parms[vi,,] <- transOrigPopt.default( normpopt[vi,,], poptDistr[vi] )
			if( 0 < length(vtwdemc$Y) )
				res$Y[ rownames(normpopt)[vi] ,,] <-  transOrigPopt.default( res$Y[ rownames(normpopt)[vi] ,,], poptDistr[vi] )
		}
		res
	})


# link on ?twDEMC, better delete resulting twDEMC.Rd
setMethodS3("twDEMC","default", function( 
		### Initialize and call \code{\link{twDEMCBlockInt}}. 
		Zinit, ##<< initial population: a numeric array (d x M0 x Npop) or an object of class twDEMC  
		...	   ##<< further arguments to \code{\link{twDEMCBlockInt}}
	){
		stop("unknown type of initialization argument Zinit")
		##seealso<< 
		## \code{\link{twDEMCBlockInt}}
		## \code{\link{twDEMC.array}}
		## \code{\link{twDEMC.twDEMC}}
	})


.makeTableDensityTest <- function(){
	df <- c(1:10,15,20,30,50)
	tmp <- -sapply( df, function(dfi){ signif(getRLogDenQuantile(NULL, maxLogDen=0, df=dfi),3) })
	plot( tmp~df)
	ds <- data.frame(df=df,diffLogDen=tmp)
	cat(twDf2wikiTable(ds))	# package twMisc
	lm1 <- lm(tmp~df)
	abline(lm1)
	lm2 <- lm(tmp[-(1:5)]~df[-(1:5)])
	abline(lm2, col="red")		#about 0.6 per parameter 
}


