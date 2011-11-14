logDenMultiTemp <- function(
	### LogDen function for the multi temperature test case 
	theta,			##<< the parameter vector (here scalar).
	logDenAccept=numeric(0),		##<< scalar: logDen for parms from revious run for two step Metropolis decision
	metropolisStepTemp=c(yPrior=1),	##<< numeric named vector: the temperature for internal metropolis step
	... 			##<< any other arguments passed to fModel
	,thetaPrior = 0.8	##<< the prior estimate of the parameters
	,sdTheta = 0.7		##<< the standard deviation of the prior estimate
	,theta0 = 1 
	,theta1 = 7
	,offset= -500
	,maxy=12.964
){
	.yPrior <- as.vector(logDenDs1MultiTemp( theta, thetaPrior, sdTheta))
	if( !is.na((logDenXParms <- logDenAccept["yPrior"])) & (logDenXParms<0)){
		logr = (.yPrior - logDenXParms) / metropolisStepTemp["yPrior"]
		if ( is.numeric(logr) & (logr) <= log(runif(1)) ){
			#reject
			return(c(yPrior=-Inf,y=NA))
		}
	}
	.y <- as.vector(logDenDs2MultiTemp( theta, theta0=theta0, theta1=theta1, offset=offset, maxy=maxy)) 
	c(yPrior=.yPrior,y=.y)
}


logDenDs1MultiTemp <- function(
	### Prior for theta
	theta
	,thetaPrior	
	,sdTheta
){ 
	-1/2 * ( (theta-thetaPrior)/sdTheta )^2
}

logDenDs2MultiTemp <- function(
	### Misfit for theta
	theta	##<<
	,theta0
	,theta1=7
	,offset=-100
	,maxy=numeric(0)
){ 
	sd1 <- theta0/2
	m1 <- theta0
	y1 <- dnorm(theta, mean=m1, sd=sd1 )
	y1Scaled <- y1 / dnorm(m1, mean=m1, sd=sd1 )
	
	sd2 <- theta1*1/7
	m2 <- theta1
	y2 <- dnorm(theta, mean=m2, sd=sd2 )
	y2Scaled <- y2 / dnorm(m2, mean=m2, sd=sd2 )

	#y3 <- (cos((theta-theta0)*15))
	y3 <- (cos((theta-theta0)*70))
	#y3 <- (cos((theta-theta0)*150))
	
	#y <- y1Scaled*3+y2Scaled*12+y3
	y <- y1Scaled*6+y2Scaled*15+y3
	if( 0==length(maxy) )
		maxy = max(y)
	yScaled <- 300*y/maxy-300+offset
	attributes(yScaled)$maxy <- maxy
	yScaled
}


.tmp.f <- function(
	### plotting the 2d case with need for temperated prior
){
	#thetaPrior=0.8
	#sdTheta=0.25
	thetaPrior=0.5
	sdTheta=0.28
	offset <- -500
	theta0=0.9
	theta1=10
	
	x <- seq(0,theta1+1,length.out=1600)
	yPrior <- logDenDs1MultiTemp(x,thetaPrior,sdTheta)
	#yPrior <- dnorm(x, mean=thetaPrior, sd=sdTheta)
	#yPriorScaled <- yPrior/max(yPrior)
	
	#mtrace(logDenDs2MultiTemp)
	y <- logDenDs2MultiTemp(x, theta0=theta0, theta1=theta1,offset=offset)
	yConv <- yPrior+y

	plot(c(x,x),c(y,yConv),type="n", xlab="theta", ylab="LogDensity",ylim=c(-1400,0)+offset)
	lines( x, y )
	#lines(x,yPriorScaled,col="blue")
	lines(x,yPrior+offset,col="blue")
	lines(x,yConv,col="maroon")

	plot(c(x,x),c(y,yConv),type="n", xlab="theta", ylab="LogDensity",ylim=c(-400,0)+offset,xlim=c(0,3))
	lines( x, y )
	#lines(x,yPriorScaled,col="blue")
	lines(x,yPrior+offset,col="blue")
	lines(x,yConv,col="maroon")
	
	#mtrace(logDenMultiTemp)
	yConv2 <- logDenMultiTemp(x,maxy=attributes(y)$maxy)
	checkEquals(c(yPrior=yPrior,y=y),yConv2)
	
	py <- {tmp <- exp(y); tmp/max(tmp[x<2])}
	pyPrior <- {tmp <- exp(yPrior); tmp/max(tmp)}
	pyConv <- {tmp <- exp(yConv); tmp/max(tmp)}
	
	plot(x,x,ylim=c(0,1),xlim=c(0,1.5),type="n", xlab="theta", ylab="scaled p(theta)")
	lines(x,py)
	lines(x,pyPrior,col="blue")
	abline(v=thetaPrior,col="blue")
	lines(x,pyConv,col="maroon")
	abline(v=theta0,col="maroon")
}

