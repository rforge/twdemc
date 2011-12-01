dinvchisq <- function (
	### density function and for the (scaled) inverse-chi-squared distribution.
	x					##<< vector of quantiles
	, df				##<< degrees of freedom parameter, usually represented as nu
	, scale = 1/df		##<< scale parameter, usually represented as lambda.
	, log = FALSE		##<< Logical. If log=TRUE, then the logarithm of the density is returned.
){
	# adopted from copyLeft pakcage LaplaceDemon, as not available for R 2.10 which runs on the cluster
	##seealso<< \code{\link{dinvchisq}}
	x <- as.vector(x)
	df <- as.vector(df)
	scale <- as.vector(scale)
	if (any(x <= 0)) 
		stop("x must be positive.")
	if (any(df <= 0)) 
		stop("The df parameter must be positive.")
	if (any(scale <= 0)) 
		stop("The scale parameter must be positive.")
	NN <- max(length(x), length(df), length(scale))
	x <- rep(x, len = NN)
	df <- rep(df, len = NN)
	scale <- rep(scale, len = NN)
	nu <- df/2
	dens <- nu * log(nu) - log(gamma(nu)) + nu * log(scale) - 
		(nu + 1) * log(x) - (nu * scale/x)
	if (log == FALSE) 
		dens <- exp(dens)
	dens
}
attr(dinvchisq,"ex") <- function(){
	x <- dinvchisq(1,1,1)
	x <- rinvchisq(10,1)
	
	#Plot Probability Functions
	x <- seq(from=0.1, to=5, by=0.01)
	plot(x, dinvchisq(x,0.5,1), ylim=c(0,1), type="l", main="Probability Function",
		ylab="density", col="red")
	lines(x, dinvchisq(x,1,1), type="l", col="green")
	lines(x, dinvchisq(x,5,1), type="l", col="blue")
	legend(3, 0.9, expression(paste(nu==0.5, ", ", lambda==1),
			paste(nu==1, ", ", lambda==1), paste(nu==5, ", ", lambda==1)),
		lty=c(1,1,1), col=c("red","green","blue"))	
}


rinvchisq <- function (
	### random number function and for the (scaled) inverse-chi-squared distribution.
	n				## the number of observations. If length(n) > 1, then the length is taken to be the number required.
	, df
	, scale = 1/df
){
	##seealso<< \code{\link{rinvchisq}}
	df <- rep(df, len = n)
	scale <- rep(scale, len = n)
	if (any(df <= 0)) 
		stop("The df parameter must be positive.")
	if (any(scale <= 0)) 
		stop("The scale parameter must be positive.")
	x <- (df * scale)/rchisq(n, df = df)
	return(x)
}





