## Test unit 'twDEMC'
# twUtestF(twDEMCBatch)
# twUtestF(twDEMCBatch,"test.goodStart")

.setUp <- function(){
	data(twLinreg1)
	attach( twLinreg1 )
	#cat("hello world")
	#mtrace(test.saveAndRestart)
}

.tearDown <- function(){
	detach( twLinreg1 )
	#cat(".teardown called\n")
	#mtrace.off()
}

.tmp.f <- function(){
	library(snowfall)
	sfInit(parallel=TRUE, cpus=4)
}


test.goodStart <- function(){
	# nice starting values
	.nPops=2
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior = thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	#mtrace(twDEMCBatchInt)
	#undebug(twDEMCBatchInt)
	res <-  twDEMCBatch(
		Zinit=Zinit, nGen=200, nBatch=50,
		restartFilename=NULL,
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops
	)
	str(res)
	
	rescoda <- as.mcmc.list(res)
	plot(rescoda)
	(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
	(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
	
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

#twUtestF(twDEMCBatch,"test.saveAndRestart")
test.saveAndRestart <- function(){
	# testing save and restarting from that file
	.nPops=2
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,				### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior= thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	restartFilename="runittwDEMCBatch_saveAndRestart.RData"
	unlink(restartFilename)
	checkTrue( !file.exists(restartFilename) )
	
	.nPops=2
	.thin=5
	res <-  twDEMCBatch(
		Zinit=Zinit,  
		60, 50,
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		
		controlTwDEMC = list(thin=.thin),
		T0=20,
		nGenBurnin=150,
		restartFilename=restartFilename
	)
	#str(res)
	checkEquals(60/.thin+1, nrow(res$rLogLik) )
	
	checkTrue( file.exists(restartFilename) )
	if( exists("resRestart.twDEMC")) rm( resRestart.twDEMC )
	checkTrue( !exists( "resRestart.twDEMC" ) )
	load( file=restartFilename )
	checkTrue( exists( "resRestart.twDEMC" ) )
	unlink(restartFilename)
	
	#str(resRestart.twDEMC)
	#mtrace(twDEMCBatch)
	#mtrace(twDEMCBatchInt)
	rm(.nPops)	#test if argument values referring to variables have been substituted by their values
	#names(res$batchCall)
	res2 <- twDEMCBatch(resRestart.twDEMC)		#should figure out burnin and nGen by itself 
	checkEquals(60/.thin+1, nrow(res2$rLogLik) )
	res3 <- twDEMCBatch(resRestart.twDEMC, nGen=70)
	checkEquals(70/.thin+1, nrow(res3$rLogLik) )
	res4 <- twDEMCBatch(resRestart.twDEMC, nGen=40)	
	checkEquals(50/.thin+1, nrow(res4$rLogLik) )
	res5 <- twDEMCBatch(resRestart.twDEMC, nGen=300)	#exeeding burnin	
	checkEquals(300/.thin+1, nrow(res5$rLogLik) )
	
	plot(res5$temp[,1])
	
	rescoda <- as.mcmc.list(res5) 
	plot(rescoda)
	rescoda <- as.mcmc.list(res5,start=attr(res5,"batchCall")$nGenBurnin) 
	plot(rescoda)

	(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
	(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
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
#twUtestF(twDEMCBatch,"test.saveAndRestart")

#twUtestF(twDEMCBatch,"test.badStart")
test.badStart <- function(){
	#test burnin with bad starting information
	.nPops=2
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,				### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior=thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	.sdThetaBad <- sdTheta*c(1/10,10)
	.theta0Bad <- theta0*c(10,1/10)
	Zinit <- initZtwDEMCNormal( .theta0Bad, diag(.sdThetaBad^2), nChains=8, nPops=.nPops)
	Zinit[,,1:4] <- Zinit[,,1:4] * c(2,1) 
	Zinit[,,-(1:4)] <- Zinit[,,-(1:4)] * c(1,2) 
	dim(Zinit)
	
	#mtrace(twDEMCBatchInt)
	res <-  twDEMCBatch(
		Zinit=Zinit,  
		300, 50,
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		
		#probUpDirBurnin=0.8,		
		nGenBurnin=200
	)
	str(res)
	
	rescoda <- as.mcmc.list(res)
	plot(rescoda)
	rescoda <- as.mcmc.list(res,start=attr(res,"batchCall")$nGenBurnin) 
	plot(rescoda)
	suppressWarnings({
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
	# same as goodStartSeq, but with executing debugSequential=FALSE, i.e. parallel
	.nPops=2
	argsFLogLik <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior=thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik))
	
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=8, nPops=.nPops)
	dim(Zinit)
	
	res <-  twDEMCBatch(
		200, 50,
		Zinit=Zinit,  
		fLogLik=logLikGaussian, argsFLogLik=argsFLogLik,
		nPops=.nPops,
		
		nGenBurnin=150,
		probUpDirBurnin=0.8
	)
	str(res)
	
	rescoda <- as.mcmc.list(res) 
	plot(rescoda)
	rescoda <- as.mcmc.list(res,start=attr(res,"batchCall")$nGenBurnin) 
	plot(rescoda)
	(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
	(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
	
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


