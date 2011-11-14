#------ developed in scratch postnorm2.R

normConstLaplace <- function(
	### Laplace Approximation of the normalizing constant 
	sample		##<< MCMC sample of parameters, each row of the matrix is a parameter vector
	, fLogPost 	##<< function to evaluate log of unnormalized posterior density
	, alpha=0.05	##<< quantile defining the multivariate normal ellipsoid for volume correction
	, nd = cov.rob(sample)	##<< normal density approximation: list with entries center and cov
){
	#hmu <- fPost(nd$center)
	#phimu <- dmvnorm(nd$center, mean=nd$center, sigma=nd$cov, log=FALSE)	
	#CL <- hmu/phimu
	logPhimu <- dmvnorm(nd$center, mean=nd$center, sigma=nd$cov, log=TRUE)	# normal approximation at mode
	logHmu <- fLogPost(nd$center)	# unnormalized posterior at center \hat{\theta}
	logCL <- logHmu - logPhimu
	# P: proportion of sample values inside normal ellipsoid defined by alpha
	dist2 <- mahalanobis(sample, center=nd$center, cov=nd$cov)	# squared mahalanobis each row is a vector
	distCrit2 <- qchisq(alpha, df = ncol(sample))
	P <- sum(dist2 < distCrit2)/length(dist2)
	#CLV <- CL * alpha/P
	logCLV <- logCL + log(alpha/P)
	##value<< a list with estimates of the logs of the normalizing constant of the posterior
	c(logCL=logCL		##<< Laplace estimate
		, logCLV=logCLV)	##<< volume corrected Laplace estimate
	##end<<
}
#mtrace(normConstLaplace)
attr(normConstLaplace,"ex") <- function(){
	mu = c(0,0)
	sd = c(1,2)
	corr = diag(nrow=2)
	corr[1,2] <- corr[2,1] <- 0.49
	Sigma = diag(sd, nrow=length(sd)) %*% corr %*% diag(sd,nrow=length(sd)) 
	n <- 1000
	cTrue <- 10
	sample <- rmvnorm(n,mean=mu,sigma=Sigma)
	fPost <- function(theta){ cTrue* dmvnorm(theta,mean=mu,sigma=Sigma, log=FALSE)}
	fLogPost <- function(theta){ log(cTrue) + dmvnorm(theta,mean=mu,sigma=Sigma, log=TRUE)}
	exp(normConstLaplace( sample, fLogPost ))
	
	#normd <- cov.rob(sample)
	#(res <-  boot( sample, function(sample, i){normConstLaplace(sample[i,], fLogDen=fLogDen, nd=normd)}, 100))
	#(res2 <-  boot( sample, function(sample, i){normConstLaplace(sample[i,], fLogDen=fLogDen)}, 100))
}

normConstLaplaceBridge <- function(
	### Laplace bridge Approximation of the normalizing constant
	sample		##<< MCMC sample of parameters, each row of the matrix is a parameter vector
	, logPostSample	##<< the unnormalized posteror density of the sample
	, fLogPost 		##<< function to evaluate unnormalized Log-posterior for a given parameter vector
	, alpha=0.05	##<< quantile defining the multivariate normal ellipsoid for volume correction
	, nd = cov.rob(sample)	##<< normal density approximation: list with entries center and cov
	, M = 1e4		##<< number of additional simulations of posterior
	, logCest = normConstLaplace(sample, fLogPost, alpha, nd)["logCLV"]	##<< initial estimate of the normalizing constant
){
	m <- nrow(sample)
	# draw a second sample from normal approximation
	sample2 <- rmvnorm(M,mean=nd$center,sigma=nd$cov)
	# calcalate unnormalized posterior Density for that
	logPost2 <- fLogPost(sample2)
	tmp.f <- function(){
		postSample <- exp(logPostSample)
		post2 <- exp(logPost2)
		Cest <- exp(logCest)
		gamma2 <- 1/(m*post2/Cest + M*dmvnorm(sample2,mean=nd$center,sigma=nd$cov))
		qs1 <- dmvnorm(sample,mean=nd$center,sigma=nd$cov)
		gamma1 <- 1/(m*postSample/Cest + M*dmvnorm(sample,mean=nd$center,sigma=nd$cov))
		(CMW <- sum(post2*gamma2)/M / ( sum(qs1*gamma1)/m  ))
	}
	gamma2 <- 1/(m*exp(logPost2-logCest) + M*dmvnorm(sample2,mean=nd$center,sigma=nd$cov))
	logQs1 <- dmvnorm(sample,mean=nd$center,sigma=nd$cov, log=TRUE)
	gamma1 <- 1/(m*exp(logPostSample-logCest) + M*dmvnorm(sample,mean=nd$center,sigma=nd$cov))
	# factor our maximum posterior density to avoid very small exponentials
	maxLogPost2 <- max(logPost2)
	maxLogQs1 <- max(logQs1)
	logCMW <- -maxLogPost2 + log( sum(exp(logPost2+maxLogPost2)*gamma2)/M ) - (-maxLogQs1 +log( sum(exp(logQs1+maxLogQs1)*gamma1)/m )) 
	### log of the normalization constant
}
attr(normConstLaplaceBridge,"ex") <- function(){
	mu = c(0,0)
	sd = c(1,2)
	corr = diag(nrow=2)
	corr[1,2] <- corr[2,1] <- 0.49
	Sigma = diag(sd, nrow=length(sd)) %*% corr %*% diag(sd,nrow=length(sd)) 
	n <- 1000
	cTrue <- 10
	sample <- rmvnorm(n,mean=mu,sigma=Sigma)
	fLogPost <- function(theta){ log(cTrue) + dmvnorm(theta,mean=mu,sigma=Sigma, log=TRUE)}
	logPostSample <- fLogPost(sample)
	(CMW <- exp(normConstLaplaceBridge( sample, logPostSample, logPostSample )) )
	
	#normd <- cov.rob(sample)
	#(res <-  boot( sample, function(sample, i){normConstLaplace(sample[i,], fLogDen=fLogDen, nd=normd)}, 100))
	#(res2 <-  boot( sample, function(sample, i){normConstLaplace(sample[i,], fLogDen=fLogDen)}, 100))
	
	#normd <- cov.rob(sample)
	#(res <-  boot( sample, function(sample, i){normConstLaplaceBridge(sample[i,], fLogDen=fLogDen, logDenSample, nd=normd)}, 100))
	#(res2 <-  boot( sample, function(sample, i){normConstLaplaceBridge(sample[i,], fLogDen=fLogDen, logDenSample)}, 100))
}


