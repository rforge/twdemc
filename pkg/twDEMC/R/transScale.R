.tmp.f <- function(){
	twUtestF(transOrigPopt)
}

twVarDistrVec <- function(
	### create a factor level with given variable names
	varNames	##<< character vector
){
	structure( rep("norm", length(varNames)), names=varNames )
	### a named factorvector initialized by level "norm"
}

R.methodsS3::setMethodS3("transOrigPopt","default", function( 
		### Transform vectors from normal to original scale.
	normpopt	##<< numerical vector/array with values at transformed, i.e. normal, scale
	,poptDistr=parDistr[names(normpopt),"trans"]
		### character vector/array of kind of distribution/transformation (e.g. "norm","lognorm","logitnorm")
		### values with other characters indicate no transformation
		### positions must match the positions in normpopt
	,parDistr 
		### alternative way of specifying poptDistr: 
		### dataframe with parameter names in column names and column trans, such as provided by \code{\link{twQuantiles2Coef}} 
	,...
){
	# transOrigPopt.default
	
	##seealso<<
	## \code{\link{twDEMCBlockInt}}
	## \code{\link{transNormPopt.default}}
    
    ##details<< 
    ## Often it is practical to work approximately normally distributed variables and specify mu and sigma parameters. 
    ## By simple transformations this applied for other distributions as well.
    ## For variables strictly greater than zero, the logNormal distribution can be applied.
    ## For variables bounded between zero and 1 the logitNormal distribution can be applied.
    ## 
    ## The values in \code{poptDistr} that are currently supported are: "norm", "lognorm", and "logitnorm"
	
	##details<< 
	## This generic method is provided for in several forms for first argument. \itemize{
	## \item{ as a named vector: this method  } 
	## \item{ as a matrix: \code{\link{transOrigPopt.matrix}}  } 
	## \item{ as a coda's mcmc.list: \code{\link{transOrigPopt.mcmc.list}}  }
	## \item{ as a list of class twDEMC \code{\link{transOrigPopt.twDEMC}}  }
	##}

	##details<< 
	## There are further methods deal with parameter transformations. \itemize{
	## \item{ transforming from original to normal scale: \code{\link{transNormPopt.default}}  } 
	## \item{ calculating mu and sigma at normal scale from quantiles: \code{\link{twQuantiles2Coef}}  } 
	## \item{ constraing the result list of \code{\link{twQuantiles2Coef}}, and adding variance-covariance matrix: \code{\link{twConstrainPoptDistr}}  }
	## \item{ calculating the density based on the distribution arguments \code{\link{dDistr}} }
	##}
	
	##details<< 
	## Argument \code{poptDistr} should have the same dimensions as normpopt. However, it is recycled.
	## By this way it is possible to specify only one value, or vector corresponding to the rows of a matrix.
	popt <- normpopt	#default no transformation  already normal
	bo <- (!is.na(poptDistr) & poptDistr == "lognorm"); popt[bo] <- exp(normpopt[bo])
	bo <- (!is.na(poptDistr) & poptDistr == "logitnorm"); popt[bo] <- plogis(normpopt[bo])
	popt
	### Normpopt with some values transformed by exp (poptDist=="lognorm") or plogis (poptDistr=="logitnorm").
})	
#mtrace(transOrigPopt)
attr(transOrigPopt.default,"ex") <- function(){
	upperBoundProb = 0.99	# quantile of the upper boundary
	parmsBounds = list(		# mode and upper bound
		A0 = c(10,15)		
		,D0 = c(10, 100)
		,C0 = c(0.6,0.8)
	)
	varDistr <- twVarDistrVec( names(parmsBounds) )	# by default assumed normal
	varDistr["D0"] <- "lognorm"
	varDistr["C0"] <- "logitnorm"
	parDistr <- twQuantiles2Coef( parmsBounds, varDistr, upperBoundProb=upperBoundProb, useMedian=FALSE )
	poptDistr <- twConstrainPoptDistr(c("A0","C0"), parDistr)
	all.equal( upperBoundProb, pnorm(parmsBounds$A0[2], parDistr[["A0","mu"]], parDistr[["A0","sigmaDiag"]] ) )
	
	# transform entire parameter vectors between scales
	pOrig <- transOrigPopt( parDistr$mu, parDistr=parDistr )
	# note that transform of mu slighly differs from the mode for lognormal and logitnormal 
	pOrig
	# back-transform to normal scale
	#pBack <- transNormPopt( pOrig, parDistr$trans[names(pOrig)] )	
	pBack <- transNormPopt( pOrig, parDistr=parDistr )	
	all.equal( parDistr$mu, pBack )
	
	# plot quantiles for given distributions
	pGrid <- seq(0.01,0.99,length.out=31)
	plot( qnorm(pGrid, mean=parDistr["D0","mu"], sd=parDistr["D0","sigmaDiag"]) ~ pGrid)
	plot( qlnorm(pGrid, mean=parDistr["D0","mu"], sd=parDistr["D0","sigmaDiag"]) ~ pGrid); abline(h=parmsBounds[["D0"]][1], col="grey")
	
	# plot densities for D0 parameter ranges
	dGrid <- seq(3, 80, length.out=100)
	denOrig1 <- dlnorm(dGrid, mean=parDistr["D0","mu"], sd=parDistr["D0","sigmaDiag"]) 
	plot( denOrig1 ~ dGrid, type="l"); abline(v=parmsBounds[["D0"]][1], col="grey")
	
	# now plot the same using a grid on normal scale, transforming them to original scale
	dNormGrid <- seq( parDistr["D0","mu"]-2*parDistr["D0","sigmaDiag"], parDistr["D0","mu"]+2*parDistr["D0","sigmaDiag"], length.out=100)
	dOrigGrid <- transOrigPopt(dNormGrid, parDistr["D0","trans"])
	all.equal( dNormGrid, transNormPopt(dOrigGrid, parDistr["D0","trans"]) )
	denOrig2 <- dDistr( dOrigGrid, parDistr["D0","mu"],parDistr["D0","sigmaDiag"],trans=parDistr["D0","trans"] )
	points( denOrig2 ~ dOrigGrid, col="blue" )
}

R.methodsS3::setMethodS3("transNormPopt","default", function( 
		### Transform vectors from original to normal scale.
		popt 
		### numerical vector/array with values at untransformed scale
		,poptDistr=parDistr[names(normpopt),"trans"]
		### character vector/array of kind of transformation ("lognorm"/"logitnorm")
		### values with other characters indicate no transformation
		### positions must match the positions in normpopt
		,parDistr 
		### alternative way of specifying poptDistr: 
		### dataframe with parameter names in column names and column trans, such as provided by \code{\link{twQuantiles2Coef}} 
		,...
	){
		# transNormPopt.default
		
		##seealso<<   
		## \code{\link{transOrigPopt.default}}
		## \code{\link{twDEMCBlockInt}}
		
		##details<< 
		## Argument \code{poptDistr} should have the same dimensions as normpopt. However, it is recycled
		## By this way it is possible to specify only one value, or vector corresponding to the rows of a matrix.
		normpopt <- popt	#default no transformation  already normal
		bo <- (!is.na(poptDistr) & poptDistr == "lognorm"); normpopt[bo] <- log(popt[bo])
		bo <- (!is.na(poptDistr) & poptDistr == "logitnorm"); normpopt[bo] <- qlogis(popt[bo])
		normpopt
		### Argument \code{popt} with some values transformed by log (poptDist=="lognorm") or qlogis (poptDistr=="logitnorm").
	})	
#mtrace(transOrigPopt)

R.methodsS3::setMethodS3("transOrigPopt","matrix", function( 
	### Applies \code{\link{transOrigPopt.default}} to each column of \code{normopt}.
	normpopt 
		### numerical matrx with values at transformed, i.e. normal, scale
	,poptDistr=parDistr[colnames(normpopt),"trans"]
		### character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt
		### Positions must match the positions in normpopt.
		### If given a single value, it is repeated.
	,parDistr 
		### Alternative way of specifying poptDistr: 
		### dataframe with parameter names in column names and column trans, such as provided by \code{\link{twQuantiles2Coef}} 
	,...
){
	#seealso<<   
	# \code{\link{transOrigPopt.default}}
	popt <- normpopt
	##details<< 
	## either poptDistr has names for each column name
	## or poptDistr has the same length as colnames(normpopt)
	if( 0==length(colnames(normpopt)) ){
		if( length(poptDistr) == 1) poptDistr <- rep(poptDistr, ncol(normpopt) )
		if( length(poptDistr) != ncol(normpopt) )
			stop("poptDistr must be of length 1 or same length as ncol(normpopt).")
	}else{ 
		if( !is.null(names(poptDistr)) )
			if( all(colnames(normpopt) %in% names(poptDistr)) )
				poptDistr=poptDistr[ colnames(normpopt) ]
		if( length(poptDistr) != ncol(normpopt) )
			stop("names of poptDistr must have entries for each column or poptDistr has no names but is of the same length as colnames(normpopt).")
	}
	#vi <- 1
	for( vi in 1:ncol(normpopt) ){
		popt[,vi] <- transOrigPopt.default( normpopt[,vi], poptDistr[vi] )
	}
	popt
})

R.methodsS3::setMethodS3("transOrigPopt","array", function( 
		### Applies \code{\link{transOrigPopt.default}} to each column, i.e. variable, of \code{normopt}.
		normpopt 
		### numerical array with values at transformed, i.e. normal, scale
		,poptDistr=parDistr[colnames(normpopt),"trans"]
		### character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt
		,parDistr 
		### Alternative way of specifying poptDistr: 
		### dataframe with parameter names in column names and column trans, such as provided by \code{\link{twQuantiles2Coef}} 
		,...
	){
		#seealso<<   
		# \code{\link{transOrigPopt.default}}
		popt <- normpopt
		##details<< 
		## either poptDistr has names for each column name
		## or poptDistr has the same length as colnames(normpopt)
		if( !is.null(names(poptDistr)) )
			if( all(colnames(normpopt) %in% names(poptDistr)) )
				poptDistr=poptDistr[ colnames(normpopt) ]
		if( length(poptDistr) != ncol(normpopt) )
			stop("names of poptDistr must have entries for each column or poptDistr has no names but is of the same length as colnames(normpopt).")
		for( vi in seq(along.with=colnames(normpopt)) ){
			popt[,vi,] <- transOrigPopt.default( normpopt[,vi,], poptDistr[vi] )
		}
		popt
	})


#twUtest(transOrigPopt, "test.transCoda" )
R.methodsS3::setMethodS3("transOrigPopt","mcmc.list", function( 
	### Applies \code{\link{transOrigPopt.default}} to each entry of \code{normopt}.
	normpopt
	### numerical matrx with values at transformed, i.e. normal, scale
	,poptDistr=parDistr[ colnames(normpopt[[1]]$parms),"trans" ]
	### character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt
	,parDistr 
	### Alternative way of specifying poptDistr: 
	### dataframe with parameter names in column names and column trans, such as provided by \code{\link{twQuantiles2Coef}} 
	,...
){
	##seealso<<   
	## \code{\link{transOrigPopt.default}}
	mcmcListApply( normpopt, transOrigPopt.matrix, poptDistr=poptDistr)
})


R.methodsS3::setMethodS3("transOrigPopt","twDEMC", function( 
	### Applies \code{\link{transOrigPopt.default}} to each column of parameters in \code{vtwdemc}.
	normpopt
	### numerical matrx with values at transformed, i.e. normal, scale
	,poptDistr=parDistr[ colnames(normpopt$parms),"trans" ]
	### character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt
	,parDistr 
	### Alternative way of specifying poptDistr: 
	### dataframe with parameter names in column names and column trans, such as provided by \code{\link{twQuantiles2Coef}} 
	,...
){
	##seealso<<   
	## \code{\link{transOrigPopt.default}}
	res <- normpopt
	res$parms <- transOrigPopt.array(normpopt$parms, poptDistr)
	# in addition transform parameters in Y
	if( 0<length(normpopt$Y) ){
		varNames <- colnames(res$parms)
		res$Y[,varNames,] <- transOrigPopt.array(normpopt$Y[,varNames, ,drop=FALSE], poptDistr=poptDistr)
	}
	res
})

twCoefLnorm <- function( 
	### Calculates mu and sigma of the lognormal distribution from median and upper quantile.
	median	##<< geometric mu (median at the original exponential scale)
	,quant	##<< value at the upper quantile, i.e. practical maximum
	,sigmaFac=qnorm(0.99) 	##<< sigmaFac=2 is 95% sigmaFac=2.6 is 99% interval 
){
	##seealso<< 
	## \code{\link{twQuantiles2Coef}}
	## \code{\link{transOrigPopt.default}}
	
	# twCoefLognorm
	# mu_geo = exp(mu); sigma_geo = exp(sigma)
	# logSpace cf-Interval: mu +- n sigma (n=1.96 for cf95)
	# expSpace cf-Interval: mu_geo */ simga_geo^n
	# given geometric mu (median at the original exponential scale) and mu+sigmaFac*sigma
	tmp.mu = log(median)
	cbind(mu=tmp.mu, sigma=1/sigmaFac*(log(quant) - tmp.mu) )
	### named numeric vector: mu and sigma parameter of the lognormal distribution.
}

twCoefLnormCi <- function( 
	### Calculates mu and sigma of the lognormal distribution from lower and upper quantile, i.e. conidence interval.
	lower	##<< value at the lower quantile, i.e. practical minimum
	,upper	##<< value at the upper quantile, i.e. practical maximum
	,sigmaFac=qnorm(0.99) 	##<< sigmaFac=2 is 95% sigmaFac=2.6 is 99% interval
	,isTransScale = FALSE ##<< if true lower and upper are already on log scale
){	
	##seealso<< 
	## \code{\link{twQuantiles2Coef}}
	## \code{\link{transOrigPopt.default}}
	if( !isTRUE(isTransScale) ){
		lower <- log(lower)
		upper <- log(upper)
	}
	sigma <- (upper-lower)/2/sigmaFac
	cbind( mu=upper-sigmaFac*sigma, sigma=sigma )
	### named numeric vector: mu and sigma parameter of the lognormal distribution.
}
attr(twCoefLnormCi,"ex") <- function(){
	mu=2
	sd=c(1,0.8)
	p=0.99
	lower <- l <- qlnorm(1-p, mu, sd )		# p-confidence interval
	upper <- u <- qlnorm(p, mu, sd )		# p-confidence interval
	cf <- twCoefLnormCi(lower,upper)	
	all.equal( cf[,"mu"] , c(mu,mu) )
	all.equal( cf[,"sigma"] , sd )
}

twCoefLnormMLE <- function( 
	### Calculates mu and sigma of the lognormal distribution from mode and upper quantile.
	mle			##<< numeric vector: mode at the original scale
	,quant		##<< numeric vector: value at the upper quantile, i.e. practical maximum
	,sigmaFac=qnorm(0.99) 	##<< sigmaFac=2 is 95% sigmaFac=2.6 is 99% interval 
){	
	##seealso<< 
	## \code{\link{twQuantiles2Coef}}
	## \code{\link{transOrigPopt.default}}
	
	# twCoefLognormMLE
	# solution of 
	# (1) mle=exp(mu-sigma^2)
	# (2) sigma=1/simgaFac(log(q)-mu)
	m <- log(mle)
	s2 <- sigmaFac/2
	sigma <- -s2 + sqrt( s2^2 -(-log(quant)+m) ) #
	mu <- m + sigma^2
	cbind(mu=mu, sigma=sigma )
	### numeric matrix: columns mu and sigma parameter of the lognormal distribution. 
	### Rows correspond to rows of mle and quant
}
attr(twCoefLnormMLE,"ex") <- function(){
	# example 1: a distribution with mode 1 and upper bound 5
	(thetaEst <- twCoefLnormMLE(1,5))
	all.equal( mle <- exp(thetaEst[1] -thetaEst[2]^2), 1, check.attributes = FALSE)
	
	# plot the distributions
	xGrid = seq(0,8, length.out=81)[-1]
	dxEst <- dlnorm(xGrid, meanlog=thetaEst[1], sdlog=thetaEst[2])
	plot( dxEst~xGrid, type="l",xlab="x",ylab="density"); abline(v=c(1,5),col="gray")
	
	# example 2: true parameters, which should be rediscovered
	theta0 <- c(mu=1, sigma=0.4)
	mle <- exp(theta0[1] -theta0[2]^2)
	perc <- 0.975		# some upper percentile, proxy for an upper bound
	quant <- qlnorm(perc, meanlog=theta0[1], sdlog=theta0[2])
	(thetaEst <- twCoefLnormMLE(mle,quant=quant,sigmaFac=qnorm(perc)) )
	
	#plot the true and the rediscovered distributions
	xGrid = seq(0,10, length.out=81)[-1]
	dx <- dlnorm(xGrid, meanlog=theta0[1], sdlog=theta0[2])
	dxEst <- dlnorm(xGrid, meanlog=thetaEst[1], sdlog=thetaEst[2])
	plot( dx~xGrid, type="l")
	lines( dxEst ~ xGrid, col="blue")	#overplots the original, coincide
	
	# example 3: explore varying the uncertainty (the upper quantile)
	x <- seq(0.01,1.2,by=0.01)
	mle = 0.2
	dx <- sapply(mle*2:8,function(q99){
			theta = twCoefLnormMLE(mle,q99,qnorm(0.99))
			dx <- dDistr(x,theta[,"mu"],theta[,"sigma"],trans="lognorm")
		})
	matplot(x,dx,type="l")
}


.ofLnormMLE <- function(
	### Objective function used by \code{\link{coefLogitnormMLE}}. 
	mu					##<< the mu parameter to estimate
	,logMle				##<< the mode of the density distribution
	,quant				##<< quantiles (values) at for percile perc
	,perc=0.975 		##<< percentiles where quant is given
){
	#mle=exp( mu-sigma^2)
	if( mu<=logMle) return(.Machine$double.xmax)
	sigma=sqrt(mu - logMle)
	tmp.predp = plnorm(quant, meanlog=mu, sdlog=sigma )
	tmp.diff = tmp.predp - perc
	sum(tmp.diff^2)
}


.twCoefLnormMLE_depr <- function(
	### Estimating coefficients of logitnormal distribution from mode and upper quantile	
	mle						##<< the mode of the density function
	,quant					##<< the quantile values
	,perc=c(0.975)			##<< the probabilites for which the quantiles were specified
	,tol=mle/1000			##<< the desired accuracy
	, ... 
){
	##seealso<< \code{\link{logitnorm}}
	# twCoefLnormMLE
	logMle=log(mle)
	tmp <- optimize( .ofLnormMLE, interval=c(logMle,log(max(quant))), logMle=logMle, quant=quant, perc=perc, tol=tol,...)
	mu = tmp$minimum
	sigma = sqrt(mu - logMle)
	c(mu=mu,sigma=as.numeric(sigma))
	### named numeric vector with estimated parameters of the logitnormal distrubtion.
	### names: \code{c("mu","sigma")}
}
#mtrace(twCoefLnormMLE)

twQuantiles2Coef <- function(
	### Calculating coefficients of transformed normal distributions from quantiles.
	parmsBounds		##<< numeric matrix (nParm x 2), each row a numeric vector of length 2 specifying mode/median/lower quantile and upper quantile value. 
		## rownames must give the variabes.
		## Alternatively a list of parameters, each entry a numeric vector of length 2 specifying mode/median/lower quantile and upper quantile value.
	,varDistr		##<< character vector identifying the distribution, i.e. transformation to normal, for each parameter
	,upperBoundProb=0.99	##<< probability for upper quantile in parmsBounds
	,useMedian=FALSE	##<< if TRUE, the first entry of parmsBounds specifies the median, insted of the mode 
	,useCi=FALSE		##<< if TRUE, the first entry of parmsBounds specifies the lower quantile, insted of the mode 
){
	##seealso<<   
	## \code{\link{twCoefLnorm}}
	## \code{\link{twCoefLnormCi}}
	## \code{\link{twCoefLnormMLE}}
	## \code{\link{transOrigPopt.default}}
	## \code{\link{twDEMCBlockInt}}
	
	#check arguments
	quant <- if( is.list(parmsBounds) ){
		do.call( rbind, parmsBounds)
	} else parmsBounds
	varNames <- rownames(quant)
	varDistr <- varDistr[varNames]		#make them the same order
	n <- length(varNames)
	if( length(varDistr) != n) 
		stop("varDistr must contain a named entry for each variable in parmsBounds.")
	boNorm <- varDistr=="norm"
	boLognorm <- varDistr=="lognorm"
	boLogitnorm <- varDistr=="logitnorm"
	if( (sum(boNorm)+sum(boLognorm)+sum(boLogitnorm)) != n)
		stop("unknown distribution for some variables. Supported are norm, lognorm, and logitnorm")
	# row indices	
	i_qMedian=1L; i_qUpper=2L 
	#twEnumNames(quant,FALSE)  # gives problems with check, define variables by hand.
	upperBoundSigmaFac = qnorm(upperBoundProb)		
	mu=structure(numeric(n),names=varNames); sigmaDiag=structure(numeric(n), names=varNames)
	#
	if( any(boNorm)){
		quantTmp <- quant[boNorm, ,drop=FALSE]
		if( useCi ){
			#sigma*sigmaFac is halv the ci
			halfWidth <- (quantTmp[,i_qUpper]-quantTmp[,i_qMedian])/2
			mu[boNorm] <- quantTmp[,i_qMedian] + halfWidth 
			sigmaDiag[boNorm] <- halfWidth/upperBoundSigmaFac						
		}else{
			mu[boNorm] <- quantTmp[,i_qMedian]
			sigmaDiag[boNorm] <- (quantTmp[,i_qUpper]-quantTmp[,i_qMedian])/upperBoundSigmaFac
		}
	}
	#
	if( any(boLognorm)){
		quantTmp <- quant[boLognorm, ,drop=FALSE]
		coefs <- if( useMedian )
			twCoefLnorm(median=quantTmp[,i_qMedian], quant=quantTmp[,i_qUpper],sigmaFac=upperBoundSigmaFac)
		else if( useCi )
			twCoefLnormCi(lower=quantTmp[,i_qMedian], upper=quantTmp[,i_qUpper],sigmaFac=upperBoundSigmaFac)
		else
			twCoefLnormMLE(mle=quantTmp[,i_qMedian], quant=quantTmp[,i_qUpper],sigmaFac=upperBoundSigmaFac)
		mu[boLognorm] <- coefs[,1]
		sigmaDiag[boLognorm] <- coefs[,2]
	}
	#
	if( any(boLogitnorm) ){
		quantTmp <- quant[boLogitnorm, ,drop=FALSE]
		coefs <- if( useMedian )
				twCoefLogitnorm(median=quantTmp[,i_qMedian], quant=quantTmp[,i_qUpper],perc=upperBoundProb)
			else if (useCi)
				twCoefLogitnormCi(lower=quantTmp[,i_qMedian], upper=quantTmp[,i_qUpper],sigmaFac=upperBoundSigmaFac)
			else
				twCoefLogitnormMLE(mle=quantTmp[,i_qMedian], quant=quantTmp[,i_qUpper],perc=upperBoundProb)
		mu[boLogitnorm] <- coefs[,1]
		sigmaDiag[boLogitnorm] <- coefs[,2]
	}
	#
	##value<< parameter distribution information, dataframe with columns
	parDistr <- data.frame(
		#varName = varNames		##<< character vector: type of distribtution (norm,lognorm,logitnorm)
		trans = varDistr		##<< character vector: type of distribtution (norm,lognorm,logitnorm)
		,mu = mu				##<< numeric vector: distribution parameter mu, i.e. expected values at normal scale
		,sigmaDiag = sigmaDiag	##<< numeric vector: standard deviation for each parameter, i.e. sqrt(diagonal of matrix parameter sigma) multivariate distrubtion without correlations.
		#,invsigma = diag()		##<< numeric matrix: inverse of the distribution parameter sigma on multivariate (transformed) normal distribution
		)
		##end<<
		##<< and rownames corresponding to the parameters
	rownames(parDistr) <- varNames
	parDistr
}

twConstrainPoptDistr <- function(
	### Constrain the information on parameters to selected parameters and add variance-covariance matrix and its inverse.
	parNames		##<< subsets of parameters: either character string, or indices, or boolean vector
	,parDistr		##<< dataframe with columns trans, mu, and sigmaDiag, and variable names in rownames as returned by \code{link{twQuantiles2Coef}}
	,corrMat=NULL	##<< correlations between parameters
){
	if( !all(is.finite(match( parNames, rownames(parDistr) ))))
		warning("twConstrainPoptDistr: not all parNames in rownames(parDistr).")
	if( is.matrix(corrMat))
		stop("correlations not implemented yet.")
	#colName <- colnames(parDistr)[2]
	poptDistr <- lapply( colnames(parDistr), function(colName){ structure(parDistr[parNames,colName], names=parNames) })
	names(poptDistr) <- colnames(parDistr)
	var <- poptDistr$sigmaDiag^2
	poptDistr$sigma <- diag(var, nrow=length(var))
	poptDistr$invSigma <- diag(1/var, nrow=length(var))
	rownames(poptDistr$sigma) <- colnames(poptDistr$sigma) <- rownames(poptDistr$invSigma) <- colnames(poptDistr$invSigma) <- parNames
	### list with entries trans, mu, sigmaDiag, sigma, and sigmaInv 
	poptDistr
}

dDistr <- function(
	### Calculate density for transform of normal distribution.
	x			##<< numeric vector of quantile at original scale for to calculate density		
	,mu			##<< numeric vector (recycled)			
	,sigma		##<< numeric vector (recycled)
	,trans		##<< factor  vector: the Transformation to use levels (norm,lognorm,logitnorm)
){
	##seealso<<   
	## \code{\link{transOrigPopt.default}}
	## \code{\link{twDEMCBlockInt}}
	
	##details<< 
	## To evaluate density at original, i.e. lognormal, logitnormal, scale
	## the density at transformed normal scale has to be multiplied with the 
	## Jacobian, i.e the derivative of the transformation
	if( !is.factor(trans) ) trans <- factor(trans,levels=c("norm","lognorm","logitnorm"))
	grid <- cbind(x,mu,sigma,trans)	#recycle element
	res <- vector("numeric", length(x))
	i <- which(grid[,4]==which(levels(trans)=="norm"))
	xms <- grid[i,] 
    res[i] <- dnorm(xms[,1],mean=xms[,2],sd=xms[,3])
	i <- which(grid[,4]==which(levels(trans)=="lognorm"))
	xms <- grid[i,] 
	res[i] <- dlnorm(xms[,1],meanlog=xms[,2],sdlog=xms[,3])
	i <- which(grid[,4]==which(levels(trans)=="logitnorm"))
	xms <- grid[i,] 
	res[i] <- dlogitnorm(xms[,1],mu=xms[,2],sigma=xms[,3])
	res
	### numeric vector of length of maximum length of the arguments
}
attr(dDistr,"ex") <- function(){
	x <- seq(0.01,3,by=0.05)
	mu=0
	sigma=1
	dx <- dDistr(x,mu,sigma,trans="lognorm")
	plot( dx ~ x)
}



