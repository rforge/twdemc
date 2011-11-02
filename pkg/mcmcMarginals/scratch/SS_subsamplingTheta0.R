steppingStoneStep <- function(
	### calculate ratio of normalizing constants for sampling from power density with beta and betaPrev
	fLogDen 			##<< function of log of unnormalized density, first argument parameter vector
	,argsFLogDen=list()	##<< further arguments to fLogDen
	,fRefDen			##<< function of reference probability density
	,argsFRefDen=list()	##<< further arguments to fRefDen
	,fRefRan			##<< function to generate random draws from the reference distribution
	,argsFRefRan=list()	##<< further arguments to fRefRan
	,theta0				##<< numeric matrix: initial large sample with cases in rows
	,logDen0 = apply( theta0, 1, function(theta){ do.call(fLogDen, c(list(theta), argsFLogDen)) })			
		### logDensity of the initial sample, pass for efficiency
	,refDen0 = apply( theta0, 1, function(theta){ do.call(fRefDen, c(list(theta), argsFRefDen)) })			 
	,nSample = length(refDen0) # number of sample records to use for calulcation, defaults to using all samples
	,beta=0.1
	,betaPrev=0
	,cPrev=1
){
	nSample = length(refDen0)
	thetaCur <- do.call(fRefRan, c(list(nSample),argsFRefRan) )
	logDenCur <- do.call( fLogDen, c(list(thetaCur), argsFLogDen) )
	refDenCur <- do.call( fRefDen, c(list(thetaCur), argsFRefDen) )
	pi0 <- exp(logDenCur)^betaPrev * refDenCur^(1-betaPrev)
	
	
	.tmp.f <- function(){
		# get a sample from the reference distribution by resampling with weights from the previous sample
		# does not work because its a multiplication of densities instead of sampling from p_{k-1}
		nSample0 <- length(refDen0)
		prefDen <- exp(logDen0)^betaPrev * refDen0^(1-betaPrev)
		p0 <- logDen0/sum(logDen0)
		iCur <- sample.int(n=nSample0, size=nSample, replace=TRUE, prob=prefDen/sum(prefDen)   )
		terms1 <- exp(logDen0[iCur])/refDen0[iCur]
	}
	.tmp.f <- function(){
		plot( density(theta0) )
		ot0 <- order(theta0)
		otc <- order(thetaCur)
		lines( density(thetaCur), col="grey" )
		lines( pi0[otc] ~ thetaCur[otc], col="blue" )
	}
	terms1 <- exp(logDenCur)/pi0
	phi <- max(terms1)
	rk <- 1/nSample * phi^(beta-betaPrev) * sum(  (terms1/phi)^(beta-betaPrev) )
}
attr( steppingStoneStep, "ex") <- function(){
	fLogDen <- function(theta,logc,...){ logc+dnorm(theta,log=TRUE,...) }
	c=10
	argsFLogDen=list(logc=log(c))
	
	xGrid <- seq(-3,3,length.out=41)
	y <- sapply( xGrid, fLogDen, log(c) )
	plot( exp(y)/c ~ xGrid )
	
	if( FALSE ){ #if( require(np) ){
		bw <- npregbw(formula=y~xGrid, tol=.1, ftol=.1)
		model <- npreg(bws = bw)
		fRefDen <- function(x){ predict(x, exdat=xGrid) }
		#fRerRan <- 
	}else{
		fRefDen <- function(x){ dnorm(x, sd=1.3)}
		fRefRan <- function(n){ rnorm(n, sd=1.3)}
	}
	lines( sapply(xGrid, fRefDen) ~ xGrid, col="orange")
	
	theta0 <- matrix( rnorm(500, sd=1), ncol=1 )
	logDen0 <- apply( theta0, 1, function(theta){ do.call(fLogDen, c(list(theta), argsFLogDen)) })
	refDen0 = apply( theta0, 1, function(theta){ do.call(fRefDen, c(list(theta), argsFRefDen)) })
	hist(logDen0)
	hist(refDen0)
	
	(r01 <- steppingStoneStep( fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan, theta0=theta0 ))
	(r10 <- steppingStoneStep( fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan, theta0=theta0, beta=1 ))
}

normConstSS <- function(
	...					##<< arguments to \code{\link{SSr}} excluding \code{beta, betaPrev}, and \code{cPrev}
	,nStepStones = 10
	,betaVec = seq(0,1,length.out=nStepStones+1)		##<< values of parameter beta 
		### for each stepping stone call to \code{\link{SSr}}
		### must start with 0 and end with 1
){
	cVec <- rep(1, length(betaVec) )
	#iBeta <- 2
	#iBeta <- 3
	#iBeta <- length(betaVec)
	for( iBeta in seq_along(betaVec)[-1] ){
		beta=betaVec[iBeta]
		betaPrev=betaVec[iBeta-1]
		cPrev=cVec[iBeta-1]
		#mtrace(steppingStoneStep)
		rk <- steppingStoneStep(
				#...
				fLogDen=fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, theta0=theta0, logDen0=logDen0, refDen0=refDen0		
				,beta=beta, betaPrev=betaPrev, cPrev=cPrev
		)
		cVec[iBeta] <- cPrev*rk
	}
	cVec[iBeta]
}
attr( normConstSS, "ex") <- function(){
}
	