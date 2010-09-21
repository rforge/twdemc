twLinreg1 <- within( list(),{
	fModel <- dummyTwDEMCModel
	xval <- runif(30,min=5,max=10)
	thetaTrue = c(a=10,b=5)	# the parameter vector
	sdObs = 3*xval^0.8
	obs = fModel(thetaTrue, xval) + rnorm(length(xval), sd=sdObs)		### vector of data to compare with
	invCovar = diag(1/sdObs^2,nrow = length(sdObs))		### the inverse of the Covariance of obs (its uncertainty)
	sdTheta= thetaTrue*0.05	# 5% relativeerror
	invCovarTheta = diag(1/(sdTheta)^2,nrow = length(sdTheta))
	theta0 = thetaTrue + rnorm(length(thetaTrue),sd=sdTheta)
})

with( twLinreg1,{
	plot( xval, obs )
	abline( thetaTrue )
	abline( theta0, col="gray" )
})

save(twLinreg1,file="data/twLinreg1.RData")

