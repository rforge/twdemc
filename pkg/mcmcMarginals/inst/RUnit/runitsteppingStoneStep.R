test.sameReference <- function(){
	# providing steppingStoneStep with the same reference as the density
	# then the normalizing factor should be estimated correctly
	
	fLogDen <- function(theta,logc,...){ logc+dnorm(theta,log=TRUE,...) }
	cVal=10
	argsFLogDen=list(logc=log(cVal))
	
	xGrid <- seq(-3,3,length.out=41)
	y <- sapply( xGrid, fLogDen, log(cVal) )
	plot( exp(y)/cVal ~ xGrid )
	
	fRefDen <- function(x){ dnorm(x, sd=1)}		
	fRefRan <- function(n){ rnorm(n, sd=1)}
	lines( sapply(xGrid, fRefDen) ~ xGrid, col="orange")
	
	theta0 <- matrix( rnorm(1000, sd=1), ncol=1 )
	logDen0 <- apply( theta0, 1, function(theta){ do.call(fLogDen, c(list(theta), argsFLogDen)) })
	refDen0 = apply( theta0, 1, function(theta){ do.call(fRefDen, c(list(theta), argsFRefDen)) })
	
	#mtrace(steppingStoneStep)
	(r01 <- steppingStoneStep( fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan ))
	(r10 <- steppingStoneStep( fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan, beta=1 ))
	(r10 <- steppingStoneStep( fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan, beta=1, nSample=1000 ))
	
	r10b <- sapply( 1:30, function(i){
				(r10 <- steppingStoneStep( fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan, beta=1 ))
			})
	mR10b <- mean(r10b)
	sdR10b <- sd(r10b)
	checkInterval( cVal, mR10b-1.96*sdR10b, mR10b+1.96, msg="cVal not inside 95% cf interval")
}


test.longerTails <- function(){
	# providing steppingStoneStep with the same reference as the density
	# then the normalizing factor should be estimated correctly
	
	fLogDen <- function(theta,logc,...){ logc+dnorm(theta,log=TRUE,...) }
	cVal=10
	argsFLogDen=list(logc=log(cVal))
	
	xGrid <- seq(-3,3,length.out=41)
	y <- sapply( xGrid, fLogDen, log(cVal) )
	plot( exp(y)/cVal ~ xGrid )
	
	fRefDen <- function(x){ dnorm(x, mean=1, sd=1.3)}		
	fRefRan <- function(n){ rnorm(n, mean=1, sd=1.3)}
	lines( sapply(xGrid, fRefDen) ~ xGrid, col="orange")
	
	#mtrace(steppingStoneStep)
	(r10 <- steppingStoneStep( fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan, beta=1, nSample=1000 ))
	
	r10b <- sapply( 1:30, function(i){
				(r10 <- steppingStoneStep( fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan, beta=1, nSample=1000 ))
			})
	mR10b <- mean(r10b)
	sdR10b <- sd(r10b)
	c( mR10b-1.96*sdR10b, mR10b+1.96 )
	checkInterval( cVal, mR10b-1.96*sdR10b, mR10b+1.96, msg="cVal not inside 95% cf interval")
	
}


