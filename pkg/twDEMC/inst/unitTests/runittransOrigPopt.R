.setUp <- function(){
	data(twdemcEx1)
	suppressWarnings({
		rm( parDistr )
		rm( parmsBounds )
		rm( varDistr )
	})
	.setUpDf <- within( list(),{
		parDistrNamed = list( trans=c(b="norm", a="lognorm", c="someweired"))
		parDistr = list( trans=c("lognorm","norm"))
		parmsBounds = list(
			kS = 1/100 * c(1,100)	# differnt meanings, essentially do not constrain
			,mSA = 1 *c(1,10) 		# MM constant in SOM decomposition: in the range but larger than A0, very uncertain
		)
		varDistr <- twVarDistrVec( names(parmsBounds) )
		varDistr["kS"] <- "lognorm"
	})
	attach(.setUpDf)
}
.tearDown <- function () {
	#detach(.setUpDf)
	detach()
}


test.trans <- function(){
	#mtrace(twQuantiles2Coef)
	parDistr <- twQuantiles2Coef(parmsBounds, varDistr)
	parNames <- c("kS")
	res <- NULL; res <- twConstrainPoptDistr(parNames,parDistr )
	checkTrue( all(sapply(res, length)==1))
	checkTrue( all(c(names(parDistr),"sigma","invSigma") %in% names(res)) )
	checkEqualsNumeric( parDistr[parNames,"sigmaDiag"]^2, diag(res$sigma) )
	checkEqualsNumeric( 1/parDistr[parNames,"sigmaDiag"]^2, diag(res$invSigma) )
}


test.transCoda <- function(){
	rescoda <- as.mcmc.list(twdemcEx1)
	#mtrace(transOrigPopt.mcmc.list)
	#parDistr <- c("lognorm","norm")
	.tmp <- transOrigPopt(rescoda, parDistr$trans)

	checkTrue( is(.tmp,"mcmc.list"))
	checkEquals( exp(rescoda[[1]][,"a"]), .tmp[[1]][,"a"] )
	checkEquals( (rescoda[[1]][,"b"]), .tmp[[1]][,"b"] )
	.n <- length(rescoda)
	checkEquals( .n, length(.tmp) )
	checkEquals( exp(rescoda[[.n]][,"a"]), .tmp[[.n]][,"a"] )
	checkEquals( (rescoda[[.n]][,"b"]), .tmp[[.n]][,"b"] )
	
	.tmp2 <- transOrigPopt( rescoda, parDistrNamed$trans )
	checkEquals( .tmp, .tmp2)
}

test.transMatrix <- function(){
	#mtrace(transOrigPopt.array)
	ss <- stackChains(twdemcEx1)[,-1]
	.tmp <- transOrigPopt(ss, parDistr$trans)
	
	checkTrue( is(.tmp,"array"))
	.tmpexp <- ss
	.tmpexp[,"a"] <- exp(ss[,"a"])
	checkEquals( .tmpexp, .tmp )
}

test.transArray <- function(){
	#mtrace(transOrigPopt.array)
	ss <- twdemcEx1$pops[[1]]$parms
	.tmp <- transOrigPopt(ss, parDistr$trans)
	
	checkTrue( is(.tmp,"array"))
	.tmpexp <- ss
	.tmpexp[,"a",] <- exp(ss[,"a",])
	checkEquals( .tmpexp, .tmp )
}


test.transTwDEMC <- function(){
	#mtrace(transOrigPopt.twDEMC)
	mc1 <- concatPops(twdemcEx1)
	.tmp <- transOrigPopt(mc1, parDistr$trans)
	
	checkTrue( is(.tmp,"twDEMC"))
	.tmpexp <- mc1
	.tmpexp$parms[,"a",] <- exp(mc1$parms[,"a",])
	.tmpexp$Y[,"a",] <- exp(mc1$Y[,"a",])
	checkEquals( .tmpexp, .tmp )
	
	.tmp2 <- transOrigPopt(mc1 )	#check default automatically accessing parDistr$trans
	checkEquals( .tmp, .tmp2)
}

test.twCoefLnorm <- function(){
	q <- c(2,7)
	theta <- twCoefLnorm(q[1],q[2])
	q2 <- qlnorm(c(0.5,0.99), meanlog=theta[1], sdlog=theta[2])
	checkEqualsNumeric(q,q2)
	
	#check vectorized result
	res <- twCoefLnorm(1:5,q[2])						
	checkEquals(5,nrow(res))
	checkEquals(res[2,,drop=FALSE], theta)
	res <- twCoefLnorm(q[1],q[2]+(-2:2))
	checkEquals(5,nrow(res))
	checkEquals(res[3,,drop=FALSE], theta)
	
}

test.twCoefLnormMLE <- function(){
	q <- c(2,7)
	#mtrace(twCoefLognormMLE)
	theta <- twCoefLnormMLE(q[1],q[2])
	tmp <- plnorm(q[2], meanlog=theta[1], sdlog=theta[2])
	q2 <- qlnorm(c(0.5,0.99), meanlog=theta[1], sdlog=theta[2])
	checkEqualsNumeric(q[2],q2[2])
	mle <- exp(theta[1]-theta[2]^2)
	checkEqualsNumeric(q[1],mle)
	
	xGrid = seq(0,10, length.out=81)[-1]
	dx <- dlnorm(xGrid, meanlog=theta[1], sdlog=theta[2])
	plot( dx~xGrid, type="l")
	abline(v=c(q[1],q[2]), col=c("gray"))			#mode
	
	res <- twCoefLnormMLE(1:5,q[2])						
	checkEquals(5,nrow(res))
	checkEquals(res[2,,drop=FALSE], theta)
	res <- twCoefLnormMLE(q[1],q[2]+(-2:2))
	checkEquals(5,nrow(res))
	checkEquals(res[3,,drop=FALSE], theta)
}


.tmp.f <- function(){
	s <- qnorm(0.99)
	mu=theta[1]
	sigma=theta[2]
	logq=log(q[2])
	m = mu-sigma^2	#log(mle)
	
	identical(sigma, 1/s*( logq-mu ) )
	identical(sigma, 1/s*( logq-m-sigma^2 ) )
	identical(s*sigma,  logq-m-sigma^2  )
	identical(0,  logq-m-sigma^2-s*sigma  )
	identical(0,  -logq+m+sigma^2+s*sigma  )
	identical(0,  +sigma^2-logq+m+s*sigma  )
	all.equal(0,  +sigma^2+s*sigma-logq+m  )
	all.equal(0,  +sigma^2+s*sigma-logq+m  )
	
	all.equal(sigma,  -s/2 + sqrt( s^2/4 -(-logq+m) ))
}

.test.fitLognorm <- function(){
	theta0=c(mu=1, sigma=0.4)
	#xGrid = seq(0,qlnorm(0.999, meanlog=theta0[1], sdlog=theta0[2]), length.out=81)[-1]
	mle <- exp(theta0[1]-theta0[2]^2 )
	xGrid = seq(0,10, length.out=81)[-1]
	dx <- dlnorm(xGrid, meanlog=theta0[1], sdlog=theta0[2])
	plot( dx~xGrid, type="l")
	abline(v=exp(theta0[1]),col="gray")	#median
	abline(v=mle,col="blue")			#mode
	abline(v=exp(theta0[1]+(theta0[2]^2)/2),col="green")	#E(x), i.e. mean
	i = which.max(dx)
	xGrid[i]
	#z <- rlnorm(1e4, meanlog=theta0[1], sdlog=theta0[2])
	#plot(density(z))
	
	fLogDenLNorm <- function(
		### LogDensity of Lognormal distribution
		x
		,theta0  # the priors
	){
		mu=theta0[1]
		sd=theta0[2]
		resPrior <- ifelse( x <= 0, -Inf, {
				logx <- log(x)
				-1/2*((logx-mu)^2)/sd^2 -logx
			})
		cbind(parms=as.numeric(resPrior))
		### LogDensity of observations and prior
	}

	tmpL <- fLogDenLNorm(xGrid,theta0)
	i2 <- which.max(tmpL)
	xGrid[i2]
	plot( -tmpL ~ xGrid, ylim=c(0,2))
	abline(v=theta0[1],col="gray")
	abline(v=mle,col="blue")
	
	#thetaPrior=theta0[1];covarTheta=theta0[2]
	#mtrace(initZtwDEMCNormal)
	Zinit <- exp(initZtwDEMCNormal( theta0[1], theta0[2] ))
	#mtrace(twDEMCBlockInt)
	resMC4 <- twDEMCBatch( Zinit, nGen=800, nPops=2
		, fLogDen=fLogDenLNorm, argsFLogDen=list(theta0=theta0)
	)
	#mtrace(as.mcmc.list.twDEMC)
	plotThinned(as.mcmc.list(resMC4))
	
	res <- stackChains(thin(resMC4))
	plot( dx~xGrid, type="l")
	abline(v=mle,col="gray")			#mode
	lines( density(res[,2]), col="blue")
	abline( v=res[which.max(res[,1]),2], col="lightblue")
	# two lines match sufficiently: sampled correct log-normal distribution

	# second approach sample in log space
	fLogDenLNormLog <- function(
		### LogDensity of Lognormal distribution
		logx	 # log of the parameter
		,theta0  # the priors
	){
		mu=theta0[1]
		sd=theta0[2]
		-1/2*((logx-mu)^2)/sd^2
		### LogDensity of observations and prior
	}
	
	Zinit <- initZtwDEMCNormal( theta0[1], theta0[2] )
	#mtrace(twDEMCBlockInt)
	resMC5 <- twDEMCBatch( Zinit, nGen=800, nPops=2
		, fLogDen=fLogDenLNormLog, argsFLogDen=list(theta0=theta0)
	)
	#mtrace(as.mcmc.list.twDEMC)
	plotThinned(as.mcmc.list(resMC5))
	
	resLog <- res <- stackChains(resMC5)
	res[,2] <- exp(resLog[,2])

	plot( dx~xGrid, type="l")
	abline(v=mle,col="gray")			#mode
	lines( density(res[,2]), col="blue")
	abline( v=res[which.max(res[,1]),2], col="lightblue")	# here median
	# two lines match sufficiently: sampled correct log-normal distribution
	

	#fitting the mode
	mle <- exp(theta0[1]-theta0[2]^2)
	p <- 0.99
	q <- qlnorm(u2, meanlog=theta0[1], sdlog=theta0[2])
	#mtrace(.ofLnormMLE)
	.ofLognormMLE(theta0[1],log(mle),q,p)
	#mtrace(twCoefLnormMLE)
	(thetaEst <- twCoefLnormMLE(mle,q,tol=1e-14))
	dxEst <- dlnorm(xGrid, meanlog=thetaEst[1], sdlog=thetaEst[2])
	
	plot( dx~xGrid, type="l")
	lines( dxEst ~ xGrid, col="blue")
	qlnorm(u2, meanlog=thetaEst[1], sdlog=thetaEst[2])
}
