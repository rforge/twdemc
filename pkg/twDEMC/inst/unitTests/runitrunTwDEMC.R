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


#twUtestF(twDEMCBatch,"test.saveAndRestart")
test.updatedInvocation <- function(){
	# testing first fit with internal Metroplis step, and continuing with single Metropolis step
	.nPops=2
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,				### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		thetaPrior= thetaTrue,	### the prior estimate of the parameters
		invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=4*.nPops, nPops=.nPops)
	
	.thin=5
	argsTwDEMCBatch0 <- list(
		Zinit=Zinit  
		,nGen=60, nBatch=30
		,fLogDen=logDenGaussian, argsFLogDen=argsFLogDen
		,nPops=.nPops
		,controlTwDEMC = list(thin=.thin)
		,T0=20
		,nGenBurnin=150
		#,intResCompNames ="parms"		#better provide with twRunDEMC, if specified here results from previous twDEMC are overwritten
		,doRecordProposals=TRUE
	)
	#mtrace(twRunDEMC)
	res <- NULL; res <-do.call( twRunDEMC, list(argsTwDEMCBatch=argsTwDEMCBatch0,intResCompNames ="parms") )
	expNSteps <- 60/.thin+1
	checkEquals(expNSteps, nrow(res$rLogDen) )
	#checkEquals(c(1,4*.nPops), dim(res$resFLogDenX) )
	checkEquals(c(2,4*.nPops), dim(res$resFLogDenX) )
	#matplot((res$Y["a",,,drop=TRUE]),type="l")
	
	#mtrace(twRunDEMC)
	#mtrace(twDEMCBatch)
	#mtrace(twDEMC.twDEMC)
	#mtrace(twDEMCInt)
	#tmpf <- argsTwDEMCBatch0$fLogDen; mtrace(tmpf); argsTwDEMCBatch0$fLogDen <- tmpf #check if logKikAccept is empty
	res2 <- NULL; res2 <-do.call( twRunDEMC, list(argsTwDEMCBatch=argsTwDEMCBatch0,Zinit=res,nGen=60+30, intResCompNames=character(0)) )
	expNSteps <- (60+30)/.thin+1
	checkEquals(expNSteps, nrow(res2$rLogDen) )
	#checkEquals(0, length(res2$resFLogDenX) )	# not influencee by intResCompNames any more
	checkEquals(c(2,4*.nPops), dim(res$resFLogDenX) )
	#matplot((res2$Y["a",,,drop=TRUE]),type="l")
	
	#mtrace(twRunDEMC)
	#mtrace(twDEMCBatch)
	#mtrace(twDEMCInt)
	#mtrace(twDEMC.twDEMC)
	#tmpf <- argsTwDEMCBatch0$fLogDen; mtrace(tmpf); argsTwDEMCBatch0$fLogDen <- tmpf #check if logKikAccept has parms
	#res3 <- NULL; res3 <-do.call( twRunDEMC, list(argsTwDEMCBatch=argsTwDEMCBatch0,Zinit=res,nGen=60+30) )
}
#twUtestF(twDEMCBatch,"test.saveAndRestart")
