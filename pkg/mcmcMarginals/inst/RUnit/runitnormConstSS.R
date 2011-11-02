




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
	
	#mtrace(normConstSS)
	r1 <-  normConstSS(
			fLogDen=fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan
			,nStepStones=1		# just one step
	)
	checkEquals( cVal, r1, "wrong match of normalizing constant with same reference prior")
	
	r2 <-  normConstSS(
			fLogDen=fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan
			,nStepStones=2		# two steps
	)
	checkEquals( cVal, r2, "wrong match of normalizing constant with same reference prior")
	
	r10 <-  normConstSS(
			fLogDen=fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, fRefRan=fRefRan
			,nStepStones=10		# just one step
	)
	checkEquals( cVal, r10, "wrong match of normalizing constant with same reference prior")
	
}

test.longerTailsReference <- function(){
	# providing steppingStoneStep with the a reference density with longer tails (higher sd)
	
	fLogDen <- function(theta,logc,...){ logc+dnorm(theta,log=TRUE,...) }
	c=10
	argsFLogDen=list(logc=log(c))
	
	xGrid <- seq(-3,3,length.out=41)
	y <- sapply( xGrid, fLogDen, log(c) )
	plot( exp(y)/c ~ xGrid )
	
	fRefDen <- function(x){ dnorm(x, sd=1.2)}		
	fRefDen <- function(x){ dnorm(x, sd=0.8)}		
	lines( sapply(xGrid, fRefDen) ~ xGrid, col="orange")
	
	theta0 <- matrix( rnorm(1000), ncol=1 )
	logDen0 <- apply( theta0, 1, function(theta){ do.call(fLogDen, c(list(theta), argsFLogDen)) })
	refDen0 <- apply( theta0, 1, function(theta){ do.call(fRefDen, c(list(theta), argsFRefDen)) })
	
	r1 <-  normConstSS(
			fLogDen=fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, theta0=theta0, logDen0=logDen0, refDen0=refDen0
			,nStepStones=1		# just one step
	)
	r1
	checkEquals( c, r1, "wrong match of normalizing constant with same reference prior")
	
	r2 <-  normConstSS(
			fLogDen=fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, theta0=theta0, logDen0=logDen0, refDen0=refDen0
			,nStepStones=2		# two steps
	)
	r2
	checkEquals( c, r2, "wrong match of normalizing constant with same reference prior")
	
	r20 <-  normConstSS(
			fLogDen=fLogDen, argsFLogDen=argsFLogDen, fRefDen=fRefDen, theta0=theta0, logDen0=logDen0, refDen0=refDen0
			,nStepStones=200		# many steps
	)
	r20
	checkEquals( c, r20, "wrong match of normalizing constant with same reference prior")
	
}

