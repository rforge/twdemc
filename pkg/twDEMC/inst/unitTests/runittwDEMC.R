## Test unit 'twDEMC'
# twUtestF(twDEMC)
# twUtestF(twDEMC, "test.goodStartSeqData")
# twUtestF(twDEMC, "test_int.goodStartSeqData.plot2d")
# twUtestF(twDEMC, "dump")

.setUp <- function(){
	data(twLinreg1)
	data(twdemcEx1)
	attach(twLinreg1)
	attach(twdemcEx1)
}

.tearDown <- function(){
	detach( twdemcEx1 )
	detach( twLinreg1 )
}

.tmp.f <- function(){
	mtrace(logDenGaussian)
	mtrace(twCalcLogDenPar)
	mtrace(.doDEMCStep)
	mtrace(.doDEMCSteps)
	mtrace(.generateXProb)
	mtrace(.generateXPropChains)
	mtrace(.xStepSnooker)
	mtrace(twDEMCInt)
	mtrace(.generateXPropThin)
	mtrace(.xStepParallel)
	twUtestF(twDEMC,"test.probUpDir")
	twUtestF(twDEMC,"test.ZinittwDEMC")	#with thinning interval
	twUtestF(twDEMC,"test.goodStartSeq")
	twUtestF(twDEMC,"test.goodStartPar")
	twUtestF(twDEMC,"test.doAppendPrevLogDen")
	twUtestF(twDEMC,"test.goodStartSeqData")
	twUtestF(twDEMC)
	
	twUtestF(twDEMCBatch)

	twUtestF()

	prevWd <- setwd(file.path("tests"))
	dumpFileBasename="doRUnitError"
	#options(error=quote(dump.frames(dumpFileBasename, TRUE)))
	load(paste(dumpFileBasename,".rda"))
	debugger(dumpFileBasename)
	setwd(prevWd)
	
	
	
}

test.twCalcLogDenPar <- function(){
	xProp <- stackChains(twdemcEx1$parms)
	#mtrace(twCalcLogDenPar)
	#resFLogDenXStacked <- stackChains(resFLogDenX)
	str(res <- twCalcLogDenPar(function(x){2*x},xProp,debugSequential=TRUE))
	
	# with resFLogDenX, vector result, providing only truncated set of parameters
	.resFLogDenX0 <- matrix(1:nrow(xProp), nrow=nrow(xProp), dimnames=list(steps=NULL,"a"))
	.resFLogDenX0[1:5,,drop=FALSE]
	#function must take two arguments
	#mtrace(twCalcLogDenPar)
	str(res <- twCalcLogDenPar(function(y1,resFLogDenAcc){y1["a"]=y1["a"]*resFLogDenAcc["a"]; y1},xProp,resFLogDenX=.resFLogDenX0, intResCompNames=colnames(.resFLogDenX0), debugSequential=TRUE))
	checkEquals( nrow(xProp), length(res$logDen) )
	.exp <- xProp; .exp[,"a"] = xProp[,"a"]*.resFLogDenX0[,"a"]
	checkEquals(.exp, res$resFLogDen)
}

test.goodStartSeqData <- function(){
	#mtrace(test_int.goodStartSeqData.plot2d)
	test_int.goodStartSeqData.plot2d()
}

test_int.goodStartSeqData.plot2d <- function(){
	lmDummy <- lm( obs ~ xval, weights=1/sdObs^2)		# results without priors
	(.expTheta <- coef(lmDummy))
	(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
	(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=c("a","b")) )
	
	# nice starting values
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		#thetaPrior = thetaTrue,	### the prior estimate of the parameters
		#invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	Zinit <- initZtwDEMCNormal( theta0, .expCovTheta, nChains=8, nPops=.nPops)
	#dim(Zinit)
	
	.nGen=100
	#mtrace(sfRemoteWrapper)
	#mtrace(.doDEMCStep)
	#mtrace(logDenGaussian)
	#mtrace(.doDEMCSteps )
	#mtrace(twDEMCInt)
	res <-  twDEMC( Zinit, nGen=.nGen, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops
		,controlTwDEMC=list(thin=8)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	)
	#str(res)
	checkEquals( 8, res$thin )
	checkEquals( (.nGen%/%res$thin)+1,nrow(res$rLogDen))

	rescoda <- as.mcmc.list(res) 
	plot(rescoda)
	#gelman.diag(rescoda)
	#summary(rescoda)
	
	.tmp.f <- function(){
		matplot( res$rLogDen,type="l" )
		tmp <- data.frame( a=as.numeric(res$Y["a",,]), b=as.numeric(res$Y["b",,]), rLogDen=as.numeric(res$Y["rLogDen",,]) )
		colorFun <- colorRampPalette(c("yellow","red"))
		levelplot( rLogDen ~ a + b, tmp, col.regions = colorFun(50))
		tmp2 <- data.frame( a=as.numeric(res$parms["a",,]), b=as.numeric(res$parms["b",,]), rLogDen=as.numeric(res$rLogDen) )
		levelplot( rLogDen ~ a + b, tmp2, col.regions = colorFun(50))
		#image( as.numeric(res$Y["a",,]), as.numeric(res$Y["b",,]), as.numeric(res$Y["rLogDen",,]))
	}
	
	#str(summary(rescoda))
	suppressWarnings({	#glm fit in summary
			summary(rescoda)$statistics[,"Mean"]
			summary(rescoda)$statistics[,"SD"]
			thetaTrue
			sdTheta
			(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
			(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
		})
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( .expSdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(.expTheta, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

test.badStartSeqData <- function(){
	lmDummy <- lm( obs ~ xval, weights=1/sdObs^2)		# results without priors
	(.expTheta <- coef(lmDummy))
	(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
	(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=c("a","b")) )
	
	# nice starting values
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		#thetaPrior = thetaTrue,	### the prior estimate of the parameters
		#invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	#Zinit <- initZtwDEMCNormal( theta0, .expCovTheta, nChains=8, nPops=.nPops)
	.sdThetaBad <- sdTheta*c(1/10,10)
	.theta0Bad <- theta0*c(10,1/10)
	Zinit <- initZtwDEMCNormal( .theta0Bad, diag(.sdThetaBad^2), nChains=8, nPops=.nPops)
	Zinit[,,1:4] <- Zinit[,,1:4] * c(2,1) 
	Zinit[,,-(1:4)] <- Zinit[,,-(1:4)] * c(1,2) 
	#dim(Zinit)
	
	.nGen=500
	#mtrace(sfRemoteWrapper)
	#mtrace(.doDEMCStep)
	#mtrace(.doDEMCSteps )
	#mtrace(logDenGaussian)
	resa <-  twDEMC( Zinit, nGen=.nGen, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops,
		controlTwDEMC=list(thin=1),
		debugSequential=TRUE
	)
	#str(res)
	checkEquals((.nGen%/%resa$thin)+1,nrow(resa$rLogDen))
	rescoda <- as.mcmc.list(resa) 
	plot(rescoda)
	
	
	res <- thin(resa, start=300)
	rescoda <- as.mcmc.list(res) 
	plot(rescoda)
	
	#gelman.diag(rescoda)
	#summary(rescoda)
	
	#str(summary(rescoda))
	suppressWarnings({	#glm fit in summary
			summary(rescoda)$statistics[,"Mean"]
			summary(rescoda)$statistics[,"SD"]
			thetaTrue
			sdTheta
			(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
			(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
		})
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( .expSdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(.expTheta, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

test.badStartSeqData1D <- function(){
	lmDummy <- lm( obs-theta0["a"] ~ xval-1, weights=1/sdObs^2)		# results without priors
	theta1d <- theta0["b"]
	(.expTheta <- coef(lmDummy))
	(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
	(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=c("b")) )
	
	
	# nice starting values
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		#thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta[2,2,drop=FALSE],	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
		,theta0=theta0
	)
	do.call( logDenGaussian, c(list(theta=theta1d),argsFLogDen))
	
	#Zinit <- initZtwDEMCNormal( theta0, .expCovTheta, nChains=8, nPops=.nPops)
	.sdThetaBad <- sdTheta["b"]*c(1/10)
	.theta0Bad <- theta0["b"]*c(10)
	Zinit <- initZtwDEMCNormal( .theta0Bad, diag(.sdThetaBad^2,nrow = length(.sdThetaBad)), nChains=8, nPops=.nPops)
	Zinit[,,1:4] <- Zinit[,,1:4] * c(2) 
	Zinit[,,-(1:4)] <- Zinit[,,-(1:4)] * c(0.5) 
	#dim(Zinit)
	
	.nGen=500
	#mtrace(sfRemoteWrapper)
	#mtrace(.doDEMCStep)
	#mtrace(.doDEMCSteps )
	#mtrace(logDenGaussian)
	resa <-  twDEMC( Zinit, nGen=.nGen, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops,
		controlTwDEMC=list(thin=1)
		#,debugSequential=TRUE
	)
	#str(res)
	checkEquals((.nGen%/%resa$thin)+1,nrow(resa$rLogDen))
	matplot( resa$parms[1,,],type="l")
	
	res <- thin(resa, start=300)
	matplot( res$parms[1,,],type="l")
	
	#gelman.diag(rescoda)
	#summary(rescoda)
	
	#str(summary(rescoda))
	suppressWarnings({	#glm fit in summary
			(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){mean(res$parms[1,,i])}))
			(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){sd(as.vector(res$parms[1,,i]))}))
		})
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( .expSdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(.expTheta, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

test.goodStartSeq <- function(){
	# nice starting values
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar/1e10,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	#dim(Zinit)
	
	.nGen=100
	#mtrace(sfRemoteWrapper)
	#mtrace(.doDEMCStep)
	#mtrace(.doDEMCSteps )
	#mtrace(logDenGaussian)
	#mtrace(twCalcLogDenPar)
	#mtrace(twDEMCInt)
	res <-  twDEMC( Zinit, nGen=.nGen, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops,
		controlTwDEMC=list(thin=1),
		debugSequential=TRUE
	)
	#str(res)
	checkEquals((.nGen%/%res$thin)+1,nrow(res$rLogDen))
	
	#windows(record=TRUE)
	rescoda <- as.mcmc.list(res)
	plot(rescoda)
	#gelman.diag(rescoda)
	#summary(rescoda)
	
	#str(summary(rescoda))
	suppressWarnings({	#glm fit in summary
		summary(rescoda)$statistics[,"Mean"]
		summary(rescoda)$statistics[,"SD"]
		thetaTrue
		sdTheta
		(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
		(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
	})
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( sdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(thetaTrue, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

test.goodStartPar <- function(){
	# same as goodStartSeq, but with executing debugSequential=FALSE, i.e. parallel
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	res <-  twDEMC( Zinit, nGen=100, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops
		#fLogDenScale=-1/2
		#debugSequential=TRUE
	)
	str(res)
	
	rescoda <- as.mcmc.list(res) 
	suppressWarnings({	#glm fit in summary
		(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
		(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
	})
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( sdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(thetaTrue, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

test.upperParBounds <- function(){
	# same as goodStartPar, but giving an upper parameter bound
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	asplit <- 10.8
	Zinit0 <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	bo.keep <- as.vector(Zinit0["a",,]) <= asplit
	#mtrace(replaceZinitCases)
	Zinit1 <- replaceZinitCases(Zinit0, bo.keep)
	Zinit2 <- replaceZinitCases(Zinit0, !bo.keep)
	dim(Zinit1)
	# replace all the cases with upperParBounds
	
	#mtrace(.doDEMCStep)
	res <-  res0 <- twDEMC( Zinit0, nGen=100, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops
		#,upperParBounds=c(a=asplit)
		#,debugSequential=TRUE
		#,controlTwDEMC=list( DRgamma=0.1, minPCompAcceptTempDecr=0.16) 
	)
	res <-  res1 <- twDEMC( Zinit1, nGen=64*4, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops
		,upperParBounds=rep( list(c(a=asplit)), .nPops)
		,debugSequential=TRUE
		#,controlTwDEMC=list( DRgamma=0.1, minPCompAcceptTempDecr=0.16) 
	)
	res <-  res2 <- twDEMC( Zinit2, nGen=64*4, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops
		,lowerParBounds=rep( list(c(a=asplit)), .nPops )
		,debugSequential=TRUE
	)
	
	ss0 <- stackChains(res0)
	# importance sampling due to integrated probability
	ss1 <- stackChains(res1)
	lw1 <- twLogSumExp(ss1[,1])
	ss2 <- stackChains(res2)
	lw2 <- twLogSumExp(ss2[,1])
	lwSum <- twLogSumExp( c(lw1,lw2) )
	pSubs <- c( exp(lw1-lwSum),  exp(lw2-lwSum) )
	nSubs <- floor(pSubs/max(pSubs) * nrow(ss0))
	ssc <- rbind( ss1[sample.int(nrow(ss1),nSubs[1]),]
		, ss2[sample.int(nrow(ss2),nSubs[2]),]
	)
	
	.tmp.f <- function(){
		str(res)
		
		rescoda <- as.mcmc.list(res)
		plot(rescoda, smooth=FALSE)
	
		plot( density(ssc[,"a"]))	
		lines(density(ss0[,"a"]), col="blue")
		lines(density(ss1[,"a"]), col="green")
		lines(density(ss2[,"a"]), col="darkgreen")
		
		gelman.diag(res1)
		gelman.diag(res2)
	}
	checkTrue( all( ss1[,"a"] <= asplit ) )
	checkTrue( all( ss2[,"a"] >= asplit ) )
	
	.popmean <- colMeans( ssc[,-1])
	.popsd <- apply(ssc[,-1], 2, sd )
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( sdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(thetaTrue, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}


tmp.f <- function(){
	twUtestF(twDEMC, "test.doAppendPrevLogDen")
	mtrace(twDEMCInt)
}

test.ZinittwDEMC <- function(){
	#mtrace(inner.ZinittwDEMC)
	#test.ZinittwDEMC()
	#twUtestF(twDEMC,"test.ZinittwDEMC")	#with thinning interval
	inner.ZinittwDEMC()
}



inner.ZinittwDEMC <- function(){
	# same as goodStartSeq, but with executing debugSequential=FALSE, i.e. parallel
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,				### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior=thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	.thin=5
	res <- res0 <-  twDEMC( Zinit, nGen=50, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops,
		controlTwDEMC = list(thin=.thin),
		debugSequential=TRUE
	)
	str(res0)
	checkEquals(50/.thin+1, nrow(res0$rLogDen) )
	#mtrace(twDEMCInt)
	rm(res)
	res <- twDEMC( res0, nGen=50,
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops,
		controlTwDEMC = list(thin=.thin)
	)
	str(res)
	checkEquals(100/.thin+1, nrow(res$rLogDen) )
	
	rescoda <- as.mcmc.list(res)
	plot(rescoda)
	suppressWarnings({	#glm fit in summary
		(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
		(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
	})
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( sdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(thetaTrue, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

test.probUpDir <- function(){
	#mtrace(inner.probUpDir)
	innertest.probUpDir()
	#twUtestF(twDEMC,"test.probUpDir")
}
	
innertest.probUpDir <- function(){
	# same as goodStartSeq, but with executing debugSequential=FALSE, i.e. parallel
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	res <-  twDEMC( Zinit, nGen=100, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops,
		#fLogDenScale=-1/2,
		controlTwDEMC = list(probUpDir=0.8)
	)
	str(res)
	
	rescoda <- as.mcmc.list(res) 
	suppressWarnings({	#glm fit in summary
	(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
	(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
	})
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( sdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(thetaTrue, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

test.ofMulti <- function(){
	# same as goodStartSeq, but with executing debugSequential=FALSE, i.e. parallel
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	res <-  twDEMC( Zinit, nGen=100, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops
		#fLogDenScale=-1/2
		#debugSequential=TRUE
	)
	str(res)
	
	rescoda <- as.mcmc.list(res) 
	
	suppressWarnings({	#glm fit in summary
	(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
	(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
	})	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( sdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(thetaTrue, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}


test.doAppendPrevLogDen <- function(){
	# same as goodStartSeq, but with executing debugSequential=FALSE, i.e. parallel
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)

	res <-  twDEMC( Zinit, nGen=100, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops,
		#fLogDenScale=-1/2,
		debugSequential=TRUE
	)
	str(res)
	
	rescoda <- as.mcmc.list(res) 
	suppressWarnings({	#glm fit in summary
	(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
	(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
	})	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( sdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(thetaTrue, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

test.dump <- function(){
	# same as goodStartSeq, but with executing debugSequential=FALSE, i.e. parallel
	.nPops=2
	argsFLogDen <- list(
		fModel=function(parms){ stop("test throwing an error")},		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta	### the inverse of the Covariance of the prior parameter estimates
		#xval=xval
	)
	#do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	.nChains=8
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=.nChains, nPops=.nPops)
	dim(Zinit)
	
	body <- expression( res <-  twDEMC( Zinit, nGen=100, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops,
		#fLogDenScale=-1/2,
		debugSequential=TRUE
	))
	checkException( eval(body))

	suppressWarnings(dir.create("tmp"))
	.remoteDumpfileBasename=file.path("tmp","testDump")
	.remoteDumpfile <- paste(.remoteDumpfileBasename,".rda",sep="")

	# dump on initial parallel calculation
	unlink(.remoteDumpfile)
	checkTrue( !file.exists(.remoteDumpfile))
	#mtrace(.doDEMCStep)
	#mtrace(twCalcLogDenPar)
	body <- expression( 
		res <-  twDEMC( Zinit, nGen=100, 
			fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
			nPops=.nPops,
			remoteDumpfileBasename=.remoteDumpfileBasename,
			#fLogDenScale=-1/2,
			debugSequential=TRUE
		)
	)
	checkException( eval(body))
	checkTrue( file.exists(.remoteDumpfile))

	#same in parallel mode
	unlink(.remoteDumpfile)
	checkTrue( !file.exists(.remoteDumpfile))
	#mtrace(.doDEMCStep)
	#mtrace(twCalcLogDenPar)
	body <- expression( 
		res <-  twDEMC( Zinit, nGen=100, 
			fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
			nPops=.nPops,
			remoteDumpfileBasename=.remoteDumpfileBasename
			#fLogDenScale=-1/2,
			#debugSequential=TRUE
		)
	)
	checkException( eval(body))
	checkTrue( file.exists(.remoteDumpfile))
	
	# deprecated: dump on Metropolis step (provided initial X and logDenX)
	# because fLogDen is evaluated at least once in twDEMC before invoking .doDEMCSteps
	.tmp.f <- function(){
		unlink(.remoteDumpfile)
		checkTrue( !file.exists(.remoteDumpfile))
		#mtrace(.doDEMCStep)
		body <- expression( 
			res <-  twDEMC( Zinit, nGen=100, 
				fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
				nPops=.nPops,
				remoteDumpfileBasename=.remoteDumpfileBasename,
				#fLogDenScale=-1/2,
				debugSequential=TRUE
				,X = matrix(1,nrow=1,ncol=.nChains)
				,logDenX = 5
				,resFLogDenX = character(0)
			)
		)
		checkException( eval(body))
		checkTrue( file.exists(.remoteDumpfile))
		
		# dump on Metropolis step (provided initial X and logDenX)
		unlink(.remoteDumpfile)
		checkTrue( !file.exists(.remoteDumpfile))
		#mtrace(.doDEMCStep)
		#options(error=recover)
		body <- expression({
				#mtrace(twCalcLogDenPar);
				res <-  twDEMC( Zinit, nGen=100 
					,fLogDen=logDenGaussian, argsFLogDen=argsFLogDen
					,nPops=.nPops
					,remoteDumpfileBasename=.remoteDumpfileBasename
					#fLogDenScale=-1/2,
					#debugSequential=TRUE,
					,X = matrix(1,nrow=1,ncol=.nChains)
				)
			})
		checkException( eval(body))
		checkTrue( file.exists(.remoteDumpfile))
	}
	
	#load( .remoteDumpfile )
	#debugger(get(.remoteDumpfileBasename))
	#choose last column (18)
	#interactively execute commands from body of doDEMCStep
}


test.twoStepMetropolis <- function(){
	# nice starting values
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar/1e10,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	.nChains=8
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=.nChains, nPops=.nPops)
	#dim(Zinit)
	X <- adrop(Zinit[,dim(Zinit)[2],,drop=FALSE],2)	#last row of Zinit
	#mtrace(twCalcLogDenPar)
	#mtrace(logDenGaussian)
	XLogDen <- twCalcLogDenPar(fLogDen=logDenGaussian, xProp=t(X), resFLogDenX=numeric(0), intResCompNames="parms", argsFLogDen=argsFLogDen
		,debugSequential=TRUE)
		
	.nGen=200
	.thin=5
	#mtrace(sfRemoteWrapper)
	#mtrace(twDEMCInt)
	#mtrace(.doDEMCSteps )
	#mtrace(logDenGaussian)
	#mtrace(.doDEMCStep)
	resa <-  twDEMC( Zinit, nGen=.nGen
		,fLogDen=logDenGaussian, argsFLogDen=argsFLogDen
		,nPops=.nPops
		,X = X, logDenX=XLogDen$logDen, resFLogDenX=t(XLogDen$resFLogDen) 
		,intResCompNames="parms"
		,controlTwDEMC = list(thin=.thin)
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	)
	#str(resa)
	
	rescoda <- as.mcmc.list(resa) 
	plot(rescoda)
	
	res <- thin(resa, start=100)
	rescoda <- as.mcmc.list(res) 
	plot(rescoda)
	
	#gelman.diag(rescoda)
	#summary(rescoda)
	
	#str(summary(rescoda))
	suppressWarnings({	#glm fit in summary
			summary(rescoda)$statistics[,"Mean"]
			summary(rescoda)$statistics[,"SD"]
			thetaTrue
			sdTheta
			(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
			(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
		})
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( sdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(thetaTrue, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

test.twoStepMetropolisTemp <- function(){
	# nice starting values
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar/1e10,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	.nChains=8
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=.nChains, nPops=.nPops)
	#dim(Zinit)
	X <- adrop(Zinit[,dim(Zinit)[2],,drop=FALSE],2)	#last row of Zinit
	#mtrace(twCalcLogDenPar)
	#mtrace(logDenGaussian)
	XLogDen <- twCalcLogDenPar(fLogDen=logDenGaussian, xProp=t(X), resFLogDenX=numeric(0), intResCompNames="parms", argsFLogDen=argsFLogDen
		,debugSequential=TRUE)
	
	
	.nGen=200
	.thin=5
	#mtrace(sfRemoteWrapper)
	#mtrace(twDEMCInt)
	#mtrace(.doDEMCSteps )
	#mtrace(logDenGaussian)	#check metropolisStepTemp
	#mtrace(.doDEMCStep)
	resa <-  twDEMC( Zinit, nGen=.nGen
		,fLogDen=logDenGaussian, argsFLogDen=argsFLogDen
		,nPops=.nPops
		,X = X, logDenX=XLogDen$logDen, resFLogDenX=t(XLogDen$resFLogDen) 
		,debugSequential=TRUE
		,controlTwDEMC = list(thin=.thin, T0=10, Tend=1/10)
		,intResCompNames="parms"
	)
	#str(resa)
	plot(resa$temp[,1])
	checkTrue( all(resa$temp[1:(80%/%.thin),] > 1) )
	
	rescoda <- as.mcmc.list(resa) 
	plot(rescoda) #decreasing amplitude
	
	res <- thin(resa, start=100)
	rescoda <- as.mcmc.list(res) 
	plot(rescoda)
	
	#gelman.diag(rescoda)
	#summary(rescoda)
	
	#str(summary(rescoda))
	suppressWarnings({	#glm fit in summary
			summary(rescoda)$statistics[,"Mean"]
			summary(rescoda)$statistics[,"SD"]
			thetaTrue
			sdTheta
			(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
			(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
		})
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( sdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(thetaTrue, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

test.goodStartSeqMultiTemp <- function(){
	# same as goodStartSeq but with decreasing temperature and multiTemp
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar/1e10,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	#dim(Zinit)
	
	.nGen=100
	#mtrace(sfRemoteWrapper)
	#mtrace(.doDEMCStep)
	#mtrace(.doDEMCSteps )
	#mtrace(logDenGaussian)
	#mtrace(twCalcLogDenPar)
	#mtrace(twDEMCInt)
	res <-  twDEMC( Zinit, nGen=.nGen, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPops=.nPops,
		#controlTwDEMC=list(thin=1, T0=20, Tend=20, useMultiT=TRUE, TFix=c(parms=1)),
		controlTwDEMC=list(thin=1, T0=20, Tend=20, useMultiT=TRUE, Tprop=c(obs=0.5,parms=1), TFix=c(parms=1)),
		debugSequential=TRUE
	)
	#str(res)
	checkEquals((.nGen%/%res$thin)+1,nrow(res$rLogDen))
	
	rescoda <- as.mcmc.list(res) 
	plot(rescoda)
	#gelman.diag(rescoda)
	#summary(rescoda)
	
	#str(summary(rescoda))
	suppressWarnings({	#glm fit in summary
			summary(rescoda)$statistics[,"Mean"]
			summary(rescoda)$statistics[,"SD"]
			thetaTrue
			sdTheta
			(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
			(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
		})
	
	# 1/2 orders of magnitude around prescribed sd for theta
	.pop=1
	for( .pop in seq(along.with=.popsd) ){
		checkMagnitude( sdTheta, .popsd[[.pop]] )
	}
	
	# check that thetaTrue is in 95% interval 
	.pthetaTrue <- sapply(1:2, function(.pop){
			pnorm(thetaTrue, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
		})
	checkInterval( .pthetaTrue ) 
}

