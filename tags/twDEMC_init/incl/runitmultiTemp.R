## Test unit 'twDEMC'
# twUtestF(twDEMCBatch)
# twUtestF(twDEMCBatch,"test.goodStart")

.setUp <- function(){
	argsFLogLik <- argsFLogLik_ <- list(
		thetaPrior=c(theta=0.5)
		,sdTheta=c(theta=0.28)
		,offset = -500
		,theta0 = c(theta=0.9)
		,theta1 = c(theta=10)
		,maxy=numeric(0)
	)
	x <- seq(0,8,length.out=200)
	y <- logLikDs2MultiTemp(x, theta0=argsFLogLik_$theta0, offset=argsFLogLik_$offset)
	argsFLogLik$maxy <- argsFLogLik_$maxy <-  attributes(y)$maxy
	setupDf <- within(list(),{
			argsFLogLik <- argsFLogLik_
	})
	attach( setupDf )
	#cat("hello world")
	#mtrace(test.saveAndRestart)
}

.tearDown <- function(){
	detach( twLinreg1 )
	#cat(".teardown called\n")
	#mtrace.off()
}


test.multiTemp1<- function(){
	x <- seq(0,8,length.out=800)
	x1 <- seq(0.7,2.5,length.out=800)
	res1 <- do.call( logLikMultiTemp, c(list(theta=1),argsFLogLik))
	checkTrue( res1["y"] < 500 )
	#mtrace(logLikMultiTemp)
	res1 <- sapply( x, function(thetai){do.call( logLikMultiTemp, c(list(theta=thetai),argsFLogLik))})
	plot( x, res1["yPrior",], type="l" )
	plot( x, res1["y",], type="l" )
	lines( x, res1["y",]+res1["yPrior",], type="l", col="maroon" )
	abline( v=argsFLogLik$thetaPrior, col="blue")
	abline( v=argsFLogLik$theta0, col="black")
}

test.noTemp <- function(){
	.nPops=1
	.nChains=.nPops*4
	Zinit <- initZtwDEMCNormal( argsFLogLik$thetaPrior/5, diag((argsFLogLik$sdTheta)^2,nrow=length(argsFLogLik$sdTheta)), nChains=.nChains, nPops=.nPops)
	#plot(density(Zinit["theta",,]))
	str(Zinit)
	
	.nGen=200
	.thin=1
	.nGenBurnin=.nGen
	#mtrace(twDEMCInt)
	#mtrace(logLikMultiTemp)
	#mtrace(.doDEMCStep)
	res <-  twDEMCBatch( Zinit, nGen=.nGen
		,fLogLik=logLikMultiTemp, argsFLogLik=argsFLogLik
		,nPops=.nPops
		,controlTwDEMC=list(thin=.thin)
		,nGenBurnin=0, T0=1
		,doRecordProposals = TRUE
		,intResCompNames="yPrior"
		,debugSequential=TRUE
	)
	str(res)

	matplot(res$parms["theta",,],type="l")
	
	.start=100/.thin
	#matplot( res$temp[-(1:.start),], type="l")
	matplot( twDEMCPopMeans(res$pAccept[-(1:.start),], ncol(res$temp), 2), type="l" )
	matplot( res$rLogLik[-(1:.start),], type="l")
	plot(density(res$parms["theta",,]))
	
	rest <- thin(res,start=.start)
	#mtrace(ggplotDensity.twDEMC)
	p1 <- ggplotDensity.twDEMC(rest)
	#p1
	matplot(res$parms["theta",1:20,],type="l")
	matplot(res$rLogLik[1:20,],type="l")
	
	str(res$Y)
	rownames(res$Y)
	checkEqualsNumeric( res$Y["rLogLik",,], res$Y["y",,]+res$Y["yPrior",,] )
	plot( res$Y["theta",,], res$Y["rLogLik",,])
	plot( res$Y["theta",,], res$Y["rLogLik",,], xlim=c(0.85,0.95))
	plot( res$Y["theta",,], res$Y["y",,])
	plot( res$Y["theta",,], res$Y["yPrior",,])
	#works surprisingly well
}

test.tempSingle <- function(){
	.nPops=1
	.nChains=.nPops*4
	Zinit <- initZtwDEMCNormal( argsFLogLik$thetaPrior*2, diag((argsFLogLik$sdTheta)^2,nrow=length(argsFLogLik$sdTheta)), nChains=.nChains, nPops=.nPops)
	#plot(density(Zinit["theta",,]))
	str(Zinit)
	
	.nGen=200
	.thin=1
	#mtrace(twDEMCInt)
	#mtrace(logLikMultiTemp)
	#mtrace(.doDEMCStep)
	res <-  res1 <- twDEMCBatch( Zinit, nGen=.nGen
		,fLogLik=logLikMultiTemp, argsFLogLik=argsFLogLik
		,nPops=.nPops
		,controlTwDEMC=list(thin=.thin, useMultiT=FALSE)
		,nGenBurnin=.nGen*1000, T0=800
		,probUpDirBurnin=0.5
		,intResCompNames="yPrior"
		,debugSequential=TRUE
		,doRecordProposals = TRUE
	)
	str(res)
	res <- res1
	
	matplot(res$parms["theta",,],type="l")
	
	.start=100/.thin
	matplot( res$temp[,], type="l")
	matplot( twDEMCPopMeans(res$pAccept[-(1:.start),], ncol(res$temp), 2), type="l" )
	matplot( res$rLogLik[-(1:.start),], type="l")
	plot(density(res$parms["theta",,]))
	
	rest <- thin(res,start=.start)
	#mtrace(ggplotDensity.twDEMC)
	p1 <- ggplotDensity.twDEMC(rest)
	#p1
	matplot(res$parms["theta",1:20,],type="l")
	matplot(res$rLogLik[1:20,],type="l")
	
	str(res$Y)
	rownames(res$Y)
	checkEqualsNumeric( res$Y["rLogLik",,], res$Y["y",,]+res$Y["yPrior",,] )
	plot( res$Y["theta",,], res$Y["rLogLik",,])
	plot( res$Y["theta",,], res$Y["rLogLik",,], xlim=c(0.85,0.95))
	plot( res$Y["theta",,], res$Y["y",,])
	plot( res$Y["theta",,], res$Y["yPrior",,])
	#works surprisingly well

	res <- res2 <- twDEMCBatch( res1
		, nGen=.nGen+1200
		,nGenBurnin=.nGen+1000
	)
	str(res)
	
	matplot(res$parms["theta",,],type="l")
	
	i <- min(which( apply( res$temp, 1 , max ) == 1 ))
	rest <- thin(res,start=i)
	plot(density(rest$parms["theta",,]))
	plot(density(rest$rLogLik))
	

}

test.tempMult <- function(){
	.nPops=1
	.nChains=.nPops*4
	Zinit <- initZtwDEMCNormal( argsFLogLik$thetaPrior*2, diag((argsFLogLik$sdTheta)^2,nrow=length(argsFLogLik$sdTheta)), nChains=.nChains, nPops=.nPops)
	#plot(density(Zinit["theta",,]))
	str(Zinit)
	
	.nGen=200
	.thin=1
	#mtrace(twDEMCInt)
	#mtrace(logLikMultiTemp)
	#mtrace(.doDEMCStep)
	res <-  res1 <- twDEMCBatch( Zinit, nGen=.nGen
		,fLogLik=logLikMultiTemp, argsFLogLik=argsFLogLik
		,nPops=.nPops
		,controlTwDEMC=list(thin=.thin, useMultiT=TRUE)
		,nGenBurnin=.nGen*1000, T0=200
		,doRecordProposals = TRUE
		,probUpDirBurnin=0.5
		,debugSequential=TRUE
	)
	str(res)
	res <- res1
	
	matplot(res$parms["theta",,],type="l")
	
	.start=100/.thin
	matplot( res$temp[,], type="l")
	matplot( twDEMCPopMeans(res$pAccept[-(1:.start),], ncol(res$temp), 2), type="l" )
	matplot( res$rLogLik[-(1:.start),], type="l")
	plot(density(res$parms["theta",,]))
	
	str(res$Y)
	rownames(res$Y)
	checkEqualsNumeric( res$Y["rLogLik",,], res$Y["y",,]+res$Y["yPrior",,] )
	plot( res$Y["theta",,], res$Y["rLogLik",,])
	plot( res$Y["theta",,], res$Y["rLogLik",,], xlim=c(0.85,0.95))
	plot( res$Y["theta",,], res$Y["y",,])
	plot( res$Y["theta",,], res$Y["yPrior",,])
	#works surprisingly well
	
	res <- res2 <- twDEMCBatch( res1
		, nGen=.nGen+600
		,nGenBurnin=.nGen+500
	)
	str(res)
	#res <- twDEMCBatch(res, nGen=1000)
	
	i <- min(which( apply( res$temp, 1 , max ) == 1 ))
	rest <- thin(res,start=i)
	plot(density(rest$parms["theta",,]))
	plot(density(rest$rLogLik))
	
	
}


