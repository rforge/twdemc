# model example of a 2d logPosterior function

mod2dCorExLog <- function(
	### example 2d model function of normal density + sin
	theta	##<< parameter vector with names x and y
	,xval	##<< additional argument, passed by ... in logLikGaussian
	,mu		##<< mean of the normal dist 
	,Sigma	##<< covariance matrix of the normal distribution
	,c=500	##<< multiplicative constant
){
	#?rmvnorm
	# dummyTwDEMCModel
	##seealso<<   
	## \code{\link{twDEMCInt}}
	## \code{\link{logLikGaussian}}
	tmp1 <- dmvnorm(theta, mean=mu, sigma=Sigma) 
	#(sin(theta[1]*5)+1)/2 * (sin(theta[2]*5)+1)/2 * max(tmp1)
	log((sin((theta[1]+theta[2])*5)+1)/2 * tmp1 + tmp1) - c
}
attr(mod2dCorExLog,"ex") <- function(){
	mu = c(0,0)
	sd = c(2,2)
	corr = diag(nrow=2)
	corr[1,2] <- corr[2,1] <- 0.9
	Sigma = diag(sd, nrow=length(sd)) %*% corr %*% diag(sd,nrow=length(sd))
	
	
	gridx <- gridy <- seq(-4,+4,length.out=61)
	gridX <- expand.grid(gridx, gridy)
	uden <- apply( gridX, 1, mod2dCorExLog, mu=mu, Sigma=Sigma ) 
	image( gridx, gridy,  exp(matrix(uden,nrow=length(gridx))), col = rev(heat.colors(100)), xlab="x", ylab="y" )
	image( gridx, gridy,  matrix(uden,nrow=length(gridx)), col = rev(heat.colors(100)), xlab="x", ylab="y" )
	
	plot( density(uden) )
	plot( sin(gridx) ~ gridx )
}



