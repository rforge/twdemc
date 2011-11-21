.modProp <- function( theta, xval ){ theta[1]*xval }

set.seed(0815)
twTwoDenEx1 <- within( list(), {
	thetaA <- 1
	sparce <- within( list(),{
		fModel <- .modProp	# use the linear regression for both densities
		xval <- runif(5,min=7,max=10)  # only five observations
		thetaTrue <- c(a=thetaA)	# the parameter vector
		sdObs <- 0.3
		obs <- fModel(thetaTrue, xval) + rnorm(length(xval), sd=sdObs)		### vector of data to compare with
		varObs <- sdObs^2		 
		sdTheta <- thetaTrue*0.1	# 10% relative error
		varTheta <- sdTheta^2
		invCovarTheta = diag(1/(sdTheta)^2,nrow = length(sdTheta))
		theta0 <- thetaTrue + rnorm(length(thetaTrue),sd=sdTheta)
		thetaBest <- lm(obs ~ xval -1)
	})
	rich <- within( list(),{
		fModel <- dummyTwDEMCModel	# use the linear regression for both densities
		xval <- runif(500,min=7,max=10)  # only five observations
		iUpper <- sample.int(length(xval),length(xval)*0.1) 
		# in rich parameter b varies with time
		thetaTrueD <- cbind( a=thetaA, b=rep(0.9,length(xval)) ) 
		thetaTrueD[iUpper, "b"] <- 1.1 
		sdObs <- 0.3
		obs <- sapply( seq_along(xval), function(i){ 
				fModel( thetaTrueD[i,], xval[i])}) 
				+ rnorm(length(xval), sd=sdObs)		### vector of data to compare with
		varObs <- sdObs^2	
		thetaTrue <- c( a=thetaA, b=mean(thetaTrueD[,"b"]))	# effective parameter b
		sdTheta <- thetaTrue*0.1	# 10% relative error
		varTheta <- sdTheta^2
		invCovarTheta = diag(1/(sdTheta)^2,nrow = length(sdTheta))
		theta0 <- thetaTrue + rnorm(length(thetaTrue),sd=sdTheta)
		thetaBest <- lm(obs ~ xval)
	})
})

with( twTwoDenEx1$sparce,{
	plot( xval, obs )
	abline( c(0,thetaTrue) )
	abline( c(0,theta0), col="gray" )
	abline( c(0,thetaBest), col="blue" )
	c(theta0, coef(thetaBest))
})
with( twTwoDenEx1$rich,{
		plot( xval, obs )
		abline( c(thetaTrue) )
		abline( c(theta0), col="gray" )
		abline( c(thetaBest), col="blue" )
		c(theta0, coef(thetaBest))
	})


save(twTwoDenEx1,file="data/twTwoDenEx1.RData")

