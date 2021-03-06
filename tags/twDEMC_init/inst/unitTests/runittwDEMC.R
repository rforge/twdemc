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
	mtrace(logLikGaussian)
	mtrace(twCalcLogLikPar)
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
	twUtestF(twDEMC,"test.doAppendPrevLoglik")
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

test.twCalcLogLikPar <- function(){
	xProp <- stackChains(twdemcEx1$parms)
	#mtrace(twCalcLogLikPar)
	#resFLogLikXStacked <- stackChains(resFLogLikX)
	str(res <- twCalcLogLikPar(function(x){2*x},xProp,debugSequential=TRUE))
	
	# with resFLogLikX, vector result, providing only truncated set of parameters
	.resFLogLikX0 <- matrix(1:nrow(xProp), nrow=nrow(xProp), dimnames=list(steps=NULL,"a"))
	.resFLogLikX0[1:5,,drop=FALSE]
	#function must take two arguments
	#mtrace(twCalcLogLikPar)
	str(res <- twCalcLogLikPar(function(y1,resFLogLikAcc){y1["a"]=y1["a"]*resFLogLikAcc["a"]; y1},xProp,resFLogLikX=.resFLogLikX0, intResCompNames=colnames(.resFLogLikX0), debugSequential=TRUE))
	checkEquals( nrow(xProp), length(res$logLik) )
	.exp <- xProp; .exp[,"a"] = xProp[,"a"]*.resFLogLikX0[,"a"]
	checkEquals(.exp, res$resFLogLik)
}

test.goodStartSeqData <- function(){
	#mtrace(test_int.goodStartSeqData)
	test_int.goodStartSeqData.plot2d()
}

test_int.goodStartSeqData.plot2d <- function(){
	lmDummy <- lm( obs ~ xval, weights=1/sdObs^2)		# results without priors
	(.expTheta <- coef(lmDummy))
	(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
	(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=c("a","b")) )
	
	# nice starting values
	.nPops=2
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		#thetaPrior = thetaTrue,	### the prior estimate of the parameters
		#invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	Zinit <- initZtwDEMCNormal( theta0, .expCovTheta, nChains=8, nPops=.nPops)
	#dim(Zinit)
	
	.nGen=100
	#mtrace(sfRemoteWrapper)
	#mtrace(.doDEMCStep)
	#mtrace(logLikGaussian)
	#mtrace(.doDEMCSteps )
	#mtrace(twDEMCInt)
	res <-  twDEMC( Zinit, nGen=.nGen, 
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops
		,controlTwDEMC=list(thin=8)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	)
	#str(res)
	checkEquals( 8, res$thin )
	checkEquals( (.nGen%/%res$thin)+1,nrow(res$rLogLik))

	rescoda <- as.mcmc.list(res) 
	plot(rescoda)
	#gelman.diag(rescoda)
	#summary(rescoda)
	
	.tmp.f <- function(){
		matplot( res$rLogLik,type="l" )
		tmp <- data.frame( a=as.numeric(res$Y["a",,]), b=as.numeric(res$Y["b",,]), rLogLik=as.numeric(res$Y["rLogLik",,]) )
		colorFun <- colorRampPalette(c("yellow","red"))
		levelplot( rLogLik ~ a + b, tmp, col.regions = colorFun(50))
		tmp2 <- data.frame( a=as.numeric(res$parms["a",,]), b=as.numeric(res$parms["b",,]), rLogLik=as.numeric(res$rLogLik) )
		levelplot( rLogLik ~ a + b, tmp2, col.regions = colorFun(50))
		#image( as.numeric(res$Y["a",,]), as.numeric(res$Y["b",,]), as.numeric(res$Y["rLogLik",,]))
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
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		#thetaPrior = thetaTrue,	### the prior estimate of the parameters
		#invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
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
	#mtrace(logLikGaussian)
	resa <-  twDEMC( Zinit, nGen=.nGen, 
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		controlTwDEMC=list(thin=1),
		debugSequential=TRUE
	)
	#str(res)
	checkEquals((.nGen%/%resa$thin)+1,nrow(resa$rLogLik))
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
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		#thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta[2,2,drop=FALSE],	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
		,theta0=theta0
	)
	do.call( logLikGaussian, c(list(theta=theta1d),argsFLogLik))
	
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
	#mtrace(logLikGaussian)
	resa <-  twDEMC( Zinit, nGen=.nGen, 
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		controlTwDEMC=list(thin=1)
		#,debugSequential=TRUE
	)
	#str(res)
	checkEquals((.nGen%/%resa$thin)+1,nrow(resa$rLogLik))
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
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar/1e10,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	#dim(Zinit)
	
	.nGen=100
	#mtrace(sfRemoteWrapper)
	#mtrace(.doDEMCStep)
	#mtrace(.doDEMCSteps )
	#mtrace(logLikGaussian)
	#mtrace(twCalcLogLikPar)
	#mtrace(twDEMCInt)
	res <-  twDEMC( Zinit, nGen=.nGen, 
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		controlTwDEMC=list(thin=1),
		debugSequential=TRUE
	)
	#str(res)
	checkEquals((.nGen%/%res$thin)+1,nrow(res$rLogLik))
	
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
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	res <-  twDEMC( Zinit, nGen=100, 
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops
		#fLogLikScale=-1/2
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


tmp.f <- function(){
	twUtestF(twDEMC, "test.doAppendPrevLoglik")
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
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,				### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior=thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	.thin=5
	res <- res0 <-  twDEMC( Zinit, nGen=50, 
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		controlTwDEMC = list(thin=.thin),
		debugSequential=TRUE
	)
	str(res0)
	checkEquals(50/.thin+1, nrow(res0$rLogLik) )
	#mtrace(twDEMCInt)
	rm(res)
	res <- twDEMC( res0, nGen=50,
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		controlTwDEMC = list(thin=.thin)
	)
	str(res)
	checkEquals(100/.thin+1, nrow(res$rLogLik) )
	
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
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	res <-  twDEMC( Zinit, nGen=100, 
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		#fLogLikScale=-1/2,
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
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	res <-  twDEMC( Zinit, nGen=100, 
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops
		#fLogLikScale=-1/2
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


test.doAppendPrevLoglik <- function(){
	# same as goodStartSeq, but with executing debugSequential=FALSE, i.e. parallel
	.nPops=2
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)

	res <-  twDEMC( Zinit, nGen=100, 
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		#fLogLikScale=-1/2,
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
	argsFLogLik <- list(
		fModel=function(parms){ stop("test throwing an error")},		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta	### the inverse of the Covariance of the prior parameter estimates
		#xval=xval
	)
	#do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	.nChains=8
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=.nChains, nPops=.nPops)
	dim(Zinit)
	
	body <- expression( res <-  twDEMC( Zinit, nGen=100, 
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		#fLogLikScale=-1/2,
		debugSequential=TRUE
	))
	checkException( eval(body))

	.remoteDumpfileBasename="testDump"
	.remoteDumpfile <- paste(.remoteDumpfileBasename,".rda",sep="")

	# dump on initial parallel calculation
	unlink(.remoteDumpfile)
	checkTrue( !file.exists(.remoteDumpfile))
	#mtrace(.doDEMCStep)
	#mtrace(twCalcLogLikPar)
	body <- expression( 
		res <-  twDEMC( Zinit, nGen=100, 
			fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
			nPops=.nPops,
			remoteDumpfileBasename=.remoteDumpfileBasename,
			#fLogLikScale=-1/2,
			debugSequential=TRUE
		)
	)
	checkException( eval(body))
	checkTrue( file.exists(.remoteDumpfile))

	#same in parallel mode
	unlink(.remoteDumpfile)
	checkTrue( !file.exists(.remoteDumpfile))
	#mtrace(.doDEMCStep)
	#mtrace(twCalcLogLikPar)
	body <- expression( 
		res <-  twDEMC( Zinit, nGen=100, 
			fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
			nPops=.nPops,
			remoteDumpfileBasename=.remoteDumpfileBasename
			#fLogLikScale=-1/2,
			#debugSequential=TRUE
		)
	)
	checkException( eval(body))
	checkTrue( file.exists(.remoteDumpfile))
	
	# deprecated: dump on Metropolis step (provided initial X and logLikX)
	# because fLogLik is evaluated at least once in twDEMC before invoking .doDEMCSteps
	.tmp.f <- function(){
		unlink(.remoteDumpfile)
		checkTrue( !file.exists(.remoteDumpfile))
		#mtrace(.doDEMCStep)
		body <- expression( 
			res <-  twDEMC( Zinit, nGen=100, 
				fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
				nPops=.nPops,
				remoteDumpfileBasename=.remoteDumpfileBasename,
				#fLogLikScale=-1/2,
				debugSequential=TRUE
				,X = matrix(1,nrow=1,ncol=.nChains)
				,logLikX = 5
				,resFLogLikX = character(0)
			)
		)
		checkException( eval(body))
		checkTrue( file.exists(.remoteDumpfile))
		
		# dump on Metropolis step (provided initial X and logLikX)
		unlink(.remoteDumpfile)
		checkTrue( !file.exists(.remoteDumpfile))
		#mtrace(.doDEMCStep)
		#options(error=recover)
		body <- expression({
				#mtrace(twCalcLogLikPar);
				res <-  twDEMC( Zinit, nGen=100 
					,fLogLik=logLikGaussian, argsFLogLik=argsFLogLik
					,nPops=.nPops
					,remoteDumpfileBasename=.remoteDumpfileBasename
					#fLogLikScale=-1/2,
					#debugSequential=TRUE,
					,X = matrix(1,nrow=1,ncol=.nChains)
				)
			})
		checkException( eval(body))
		checkTrue( file.exists(.remoteDumpfile))
	}
	
	#load( .remoteDumpfile )
	#debugger(testDump)
	#choose last column (18)
	#interactively execute commands from body of doDEMCStep
}


test.twoStepMetropolis <- function(){
	# nice starting values
	.nPops=2
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar/1e10,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	.nChains=8
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=.nChains, nPops=.nPops)
	#dim(Zinit)
	X <- adrop(Zinit[,dim(Zinit)[2],,drop=FALSE],2)	#last row of Zinit
	#mtrace(twCalcLogLikPar)
	#mtrace(logLikGaussian)
	XLogLik <- twCalcLogLikPar(fLogLik=logLikGaussian, xProp=t(X), resFLogLikX=numeric(0), intResCompNames="parms", argsFLogLik=argsFLogLik
		,debugSequential=TRUE)
		
	.nGen=200
	.thin=5
	#mtrace(sfRemoteWrapper)
	#mtrace(twDEMCInt)
	#mtrace(.doDEMCSteps )
	#mtrace(logLikGaussian)
	#mtrace(.doDEMCStep)
	resa <-  twDEMC( Zinit, nGen=.nGen
		,fLogLik=logLikGaussian, argsFLogLik=argsFLogLik
		,nPops=.nPops
		,X = X, logLikX=XLogLik$logLik, resFLogLikX=t(XLogLik$resFLogLik) 
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
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar/1e10,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	.nChains=8
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=.nChains, nPops=.nPops)
	#dim(Zinit)
	X <- adrop(Zinit[,dim(Zinit)[2],,drop=FALSE],2)	#last row of Zinit
	#mtrace(twCalcLogLikPar)
	#mtrace(logLikGaussian)
	XLogLik <- twCalcLogLikPar(fLogLik=logLikGaussian, xProp=t(X), resFLogLikX=numeric(0), intResCompNames="parms", argsFLogLik=argsFLogLik
		,debugSequential=TRUE)
	
	
	.nGen=200
	.thin=5
	#mtrace(sfRemoteWrapper)
	#mtrace(twDEMCInt)
	#mtrace(.doDEMCSteps )
	#mtrace(logLikGaussian)	#check metropolisStepTemp
	#mtrace(.doDEMCStep)
	resa <-  twDEMC( Zinit, nGen=.nGen
		,fLogLik=logLikGaussian, argsFLogLik=argsFLogLik
		,nPops=.nPops
		,X = X, logLikX=XLogLik$logLik, resFLogLikX=t(XLogLik$resFLogLik) 
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

