.tmp.f <- function(){
	twUtestF(transOrigPopt)
}

### Distribution supported by transScale
twVarDistrLevels <- as.factor( c("norm","lognorm","logitnorm") )

twVarDistrVec <- function(
	### create a factor level with given variable names
	varNames	##<< character vector
){
	structure( rep(twVarDistrLevels[1], length(varNames)), names=varNames )
	### a named vector intiialized by factor "norm"
}

setMethodS3("transOrigPopt","default", function( 
		### Transform vectors from normal to original scale.
	normpopt, 
		### numerical vector/array with values at transformed, i.e. normal, scale
	poptDistr=eval(parse(text="parDistr$trans[names(normpopt)]")),
		### character vector/array of kind of transformation ("lognorm"/"logitnorm")
		### values with other characters indicate no transformation
		### default assumes vector parDistr$trans in environement 
	...
){
	# transOrigPopt.default
	
	##seealso<<   
	## \code{\link{twDEMCInt}}
	## \code{\link{transNormPopt.default}}
	
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
	##}
	
	##details<< 
	## Argument \code{poptDistr} should have the same dimensions as normpopt. However, it is recycled
	## By this way it is possible to specify only one value, or vector corresponding to the rows of a matrix.
	popt <- normpopt	#default no transformation  already normal
	bo <- (!is.na(poptDistr) & poptDistr == "lognorm"); popt[bo] <- exp(normpopt[bo])
	bo <- (!is.na(poptDistr) & poptDistr == "logitnorm"); popt[bo] <- plogis(normpopt[bo])
	popt
	### Normpopt with some values transformed by exp (poptDist=="lognorm") or plogis (poptDistr=="logitnorm").
})	
#mtrace(transOrigPopt)

setMethodS3("transNormPopt","default", function( 
		### Transform vectors from original to normal scale.
		popt, 
		### numerical vector/array with values at untransformed scale
		poptDistr=eval(parse(text="parDistr$trans[names(normpopt)]")),
		### character vector/array of kind of transformation ("lognorm"/"logitnorm")
		### values with other characters indicate no transformation
		### default assumes vector parDistr$trans in environement 
		...
	){
		# transNormPopt.default
		
		##seealso<<   
		## \code{\link{transOrigPopt.default}}
		## \code{\link{twDEMCInt}}
		
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

setMethodS3("transOrigPopt","matrix", function( 
	### Applies \code{\link{transOrigPopt.default}} to each column of \code{normopt}.
	normpopt, 
		### numerical matrx with values at transformed, i.e. normal, scale
	poptDistr=eval(parse(text="parDistr$trans[colnames(normpopt)]")),
		### character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt
	...
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
		popt[,vi] <- transOrigPopt.default( normpopt[,vi], poptDistr[vi] )
	}
	popt
})

setMethodS3("transOrigPopt","array", function( 
		### Applies \code{\link{transOrigPopt.default}} to each row of \code{normopt}.
		normpopt, 
		### numerical matrx with values at transformed, i.e. normal, scale
		poptDistr=eval(parse(text="parDistr$trans[rownames(normpopt)]")),
		### character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt
		...
	){
		#seealso<<   
		# \code{\link{transOrigPopt.default}}
		popt <- normpopt
		##details<< 
		## either poptDistr has names for each column name
		## or poptDistr has the same length as colnames(normpopt)
		if( !is.null(names(poptDistr)) )
			if( all(rownames(normpopt) %in% names(poptDistr)) )
				poptDistr=poptDistr[ rownames(normpopt) ]
		if( length(poptDistr) != nrow(normpopt) )
			stop("names of poptDistr must have entries for each row or poptDistr has no names but is of the same length as rownames(normpopt).")
		for( vi in seq(along.with=rownames(normpopt)) ){
			popt[vi,,] <- transOrigPopt.default( normpopt[vi,,], poptDistr[vi] )
		}
		popt
	})


#twUtest(transOrigPopt, "test.transCoda" )
setMethodS3("transOrigPopt","mcmc.list", function( 
	### Applies \code{\link{transOrigPopt.default}} to each entry of \code{normopt}.
	normpopt, 
	### numerical matrx with values at transformed, i.e. normal, scale
	poptDistr=eval(parse(text="parDistr$trans")),
	### character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt
	...
){
	##seealso<<   
	## \code{\link{transOrigPopt.default}}
	mcmcListApply( normpopt, transOrigPopt.matrix, poptDistr)
})


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
	res$parms <- transOrigPopt.array(vtwdemc$parms, poptDistr=poptDistr)
	if( 0<length(vtwdemc$Y) ){
		varNames <- rownames(res$parms)
		res$Y[varNames,,] <- transOrigPopt.array(vtwdemc$Y[varNames,, ,drop=FALSE], poptDistr=poptDistr)
	}
	res
})

twCoefLnorm <- function( 
	### Calculates mu and sigma of the lognormal distribution from median and upper quantile.
	median	##<< geometric mu (median at the original exponential scale)
	,quant	##<< value at the upper quantile, i.e. practical maximum
	,sigmaFac=qnorm(0.99) 	##<< sigmaFac=2 is 95% sigmaFac=2.6 is 99% interval 
){	
	# twCoefLognorm
	# mu_geo = exp(mu); sigma_geo = exp(sigma)
	# logSpace cf-Interval: mu +- n sigma (n=1.96 for cf95)
	# expSpace cf-Interval: mu_geo */ simga_geo^n
	# given geometric mu (median at the original exponential scale) and mu+sigmaFac*sigma
	tmp.mu = log(median)
	cbind(mu=tmp.mu, sigma=1/sigmaFac*(log(quant) - tmp.mu) )
	### named numeric vector: mu and sigma parameter of the lognormal distribution.
}

twCoefLnormMLE <- function( 
	### Calculates mu and sigma of the lognormal distribution from mode and upper quantile.
	mle			##<< numeric vector: mode at the original scale
	,quant		##<< numeric vector: value at the upper quantile, i.e. practical maximum
	,sigmaFac=qnorm(0.99) 	##<< sigmaFac=2 is 95% sigmaFac=2.6 is 99% interval 
){	
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
	tmp.predp = plnorm(quant, mean=mu, sd=sigma )
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
	parmsBounds		##<< list parameters, each entry a numeric vector of length 2 specifying mode and upper quantile value
	,varDistr		##<< character vector identifying the distribution, i.e. transformation to normal, for each parameter
	,upperBoundProb=0.99	##<< probability for upper quantile in parmsBounds
	,useMedian=FALSE	##<< if TRUE, the first entry of parmsBounds specifies the median, insted of the mode 
){
	#check arguments
	varNames <- names(parmsBounds)
	varDistr <- varDistr[varNames]		#make them the same order
	n <- length(varNames)
	if( length(varDistr) != n) 
		stop("varDistr must contain a named entry for each variable in parmsBounds.")
	boNorm <- varDistr=="norm"
	boLognorm <- varDistr=="lognorm"
	boLogitnorm <- varDistr=="logitnorm"
	if( (sum(boNorm)+sum(boLognorm)+sum(boLogitnorm)) != n)
		stop("unknown distribution for some variables. Supported are norm, lognorm, and logitnorm")
	
	tmp <- as.data.frame(parmsBounds)
	quant = cbind(qMedian = unlist(tmp[1,]), qUpper = unlist(tmp[2,])); i_qMedian=1L; i_qUpper=2L 
	#twEnumNames(quant,FALSE)  # gives problems with check, define variables by hand.
	upperBoundSigmaFac = qnorm(upperBoundProb)		
	mu=structure(numeric(n),names=varNames); sigmaDiag=structure(numeric(n), names=varNames)

	if( any(boNorm)){
		quantTmp <- quant[boNorm,,drop=FALSE]
		mu[boNorm] <- quantTmp[,i_qMedian]
		sigmaDiag[boNorm] <- (quantTmp[,i_qUpper]-quantTmp[,i_qMedian])/upperBoundSigmaFac
	}
	
	if( any(boLognorm)){
		quantTmp <- quant[boLognorm,,drop=FALSE]
		coefs <- if( useMedian )
			twCoefLnorm(median=quantTmp[,i_qMedian], quant=quantTmp[,i_qUpper],sigmaFac=upperBoundSigmaFac)
		else
			twCoefLnormMLE(mle=quantTmp[,i_qMedian], quant=quantTmp[,i_qUpper],sigmaFac=upperBoundSigmaFac)
		mu[boLognorm] <- coefs[,1]
		sigmaDiag[boLognorm] <- coefs[,2]
	}
	
	if( any(boLogitnorm) ){
		quantTmp <- quant[boLogitnorm,,drop=FALSE]
		coefs <- if( useMedian )
				twCoefLogitnorm(median=quantTmp[,i_qMedian], quant=quantTmp[,i_qUpper],perc=upperBoundProb)
			else
				twCoefLogitnormMLE(mle=quantTmp[,i_qMedian], quant=quantTmp[,i_qUpper],perc=upperBoundProb)
		mu[boLogitnorm] <- coefs[,1]
		sigmaDiag[boLogitnorm] <- coefs[,2]
	}
	
	parDistr <- list(
		trans = varDistr
		,mu = mu
		,sigmaDiag = sigmaDiag
		#,invsigma = diag()
	)
	### parameter distribution information, list entries \describe{
	### \item{trans}{character vector: type of distribtution (norm,lognorm,logitnorm)}
	### \item{mu}{numeric vector: distribution parameter mu, i.e. expected values at normal scale}
	### \item{sigmaDiag}{numeric vector: standard deviation for each parameter, i.e. diagonal of ,matrix parameter sigma multivariate distrubtion without correlations.}
	### %\item{invsigma}{numeric matrix: inverse of the distribution parameter sigma on multivariate (transformed) normal distribution}
	### }
}

twConstrainPoptDistr <- function(
	### Constrain the information on parameters to selected parameters and add variance-covariance matrix and its inverse.
	parNames		##<< subsets of parameters: either character string, or indices, or boolean vector
	,parDistr		##<< the information for all possible parameters
	,corrMat=NULL	##<< correlations between parameters
){
	#parName
	poptDistr <- list(
		trans=parDistr$trans[parNames]
		,mu= parDistr$mu[parNames]
		,sigmaDiag=parDistr$sigmaDiag[parNames]
	)
	var <- poptDistr$sigmaDiag^2
	poptDistr$sigma <- diag(var, nrow=length(var))
	poptDistr$invSigma <- diag(1/var, nrow=length(var))
	if( is.matrix(corrMat))
		stop("correlations not implemented yet.")
	poptDistr
}

dDistr <- function(
	### Calculate density for transform of normal distribution.
	x			##<< numeric vector of quantile at original scale for to calculate density		
	,mu			##<< numeric vector (recycled)			
	,sigma		##<< numeric vector (recycled)
	,trans		##<< factor  vector: the Transformation to use levels (norm,lognorm,logitnorm)
){
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



