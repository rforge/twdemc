## Test unit 'twDEMC'
# twUtestF(twDEMC)
# twUtestF(twDEMC, "test.goodStartSeqData")
# twUtestF(twDEMC, "test_int.goodStartSeqData.plot2d")
# twUtestF(twDEMC, "dump")

.setUp <- function(){
	data(twLinreg1)
	#data(twdemcEx1)
	attach(twLinreg1)
	#attach(twdemcEx1)
}

.tearDown <- function(){
	#detach( twdemcEx1 )
	detach( twLinreg1 )
}


.tmp.f <- function(){
}

test.distinctLogDen <- function(){
	lmDummy <- lm( obs ~ xval, weights=1/sdObs^2)		# results without priors
	(.expTheta <- coef(lmDummy))
	(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
	(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=c("a","b")) )
	#plot( obs ~ xval )
	#abline(lmDummy)
	#abline(thetaTrue, col="red")
	
	# nice starting values
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		#thetaPrior = thetaTrue,	### the prior estimate of the parameters
		#invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	#mtrace(logDenGaussian)
	do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	.nPop=3
	.nChainPop=4
	.thin=4
	ZinitPops <- initZtwDEMCNormal( theta0, .expCovTheta, nChainPop=4, nPop=.nPop)
	#dim(ZinitPops)
	#head(ZinitPops[,,1])
	pops0 <- list(
		pop1 <- list(
			parms = ZinitPops[1:3,,1:4,drop=FALSE]	# the first population with less initial conditions
			,nGen=8
		)
		,pop2 <- list(
			parms = ZinitPops[,,5:8,drop=FALSE]	# the first population with less initial conditions
			,nGen=100
			,T0=20
			,Tend=5
		)
		,pop3 <- list(
			parms = ZinitPops[,,9:12,drop=FALSE]	# the first population with less initial conditions
			,nGen=0
		)
	)
	#tmp <- .checkPop(pops[[1]])
	# both fLogDen compare the same model against the same observations but use different priors
	dInfoDefault <- list(
		fLogDen=logDenGaussian
		,resCompNames=c("obs","parms")
		,TFix=c(parms=1)
		,argsFLogDen = list(
			fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
			obs=obs,				### vector of data to compare with
			invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
			xval=xval
		)
	)
	dInfos0 <- list(
		logDenA = within(dInfoDefault,{
				argsFLogDen <- c( argsFLogDen, list(
						blockIndices=1
						#,thetaPrior= thetaTrue["a"]	### the prior estimate of the parameters
						# underestimates variances ,invCovarTheta = invCovarTheta[1,1,drop=FALSE]	### the inverse of the Covariance of the prior parameter estimates
						#,invCovarTheta = sdTheta^2
				)) 
			})
		,logDenB = within(dInfoDefault,{
				argsFLogDen <- c( argsFLogDen, list(
						blockIndices=2
						#,thetaPrior= thetaTrue["b"]	### the prior estimate of the parameters
						#,invCovarTheta = invCovarTheta[2,2,drop=FALSE]	### the inverse of the Covariance of the prior parameter estimates
					)) 
			})
	)
	blocks0 <- list(
		blockA <- list(compPos="a", dInfoPos="logDenA")
		,blockB <- list(compPos="b", dInfoPos="logDenB")
	)
	
	# here no nGen argument to use different nGen of pops
	#mtrace(twDEMCBlockInt)
	res <- twDEMCBlockInt( pops=pops0, dInfos=dInfos0, blocks=blocks0, controlTwDEMC=list(thin=.thin) )
	str(res$pops[[2]])
	str(res$pops[[3]])
	
	.nGenThinned=sapply(pops0, "[[", "nGen")%/%.thin+1
	for( iPop in seq_along(res$pops) ){
		resPop <- res$pops[[iPop]]
		checkEquals( .nGenThinned[iPop], nrow(resPop$parms), msg="number of states in parms mismatch" )
		checkEquals( .nGenThinned[iPop], length(resPop$temp), msg="number of states in temp mismatch" )
		checkEquals( .nGenThinned[iPop], nrow(resPop$pAccept), msg="number of states in pAccept mismatch" )
		checkEquals( .nGenThinned[iPop], nrow(resPop$resLogDen), msg="number of states in resLogDen mismatch" )
		checkEquals( .nGenThinned[iPop], nrow(resPop$logDen), msg="number of states in logDen mismatch" )
		ZinitI <- pops0[[iPop]]$parms
		checkEquals( ZinitI[nrow(ZinitI),,], resPop$parms[1,,], msg="first row of parms should correspond to the last row of parms" )		
		checkEquals( 1, resPop$temp[.nGenThinned[iPop] ], , msg="end temperature need to be 1" )		
		checkTrue( all(is.finite(resPop$parms)), msg="found non-finite values in parameters")	
		checkTrue( all(is.finite(resPop$temp)), msg="found non-finite values in temp")	
		checkTrue( all(is.finite(resPop$pAccept)), msg="found non-finite values in pAccept")	
		checkTrue( all(is.finite(resPop$resLogDen)), msg="found non-finite values in resLogDen")	
	}
	checkEquals( 1, res$pops[[1]]$temp[1], "first row of temperature of first population must be 1" )
	checkEquals( pops0[[2]]$T0, res$pops[[2]]$temp[1], "first row of temperature of second population must correspond correspond to prescribed argument" )
	
	
	.nGen=100
	#.thin=4
	#.nGen=3
	#mtrace(.checkBlock)
	#mtrace(.updateBlockTwDEMC)
	#mtrace(.updateBlocksTwDEMC)
	#mtrace(.updateIntervalTwDEMCPar)
	#mtrace(twDEMCBlockInt)
	res <- resAll <- twDEMCBlockInt( pops=pops0, dInfos=dInfos0, blocks=blocks0, nGen=.nGen, controlTwDEMC=list(thin=.thin, DRgamma=0.15)	,debugSequential=TRUE, doRecordProposals=TRUE	)
	str(res$pops[[2]]$Y)
	#windows(record=TRUE)
	plot( as.mcmc.list(res, minPopLength=10), smooth=FALSE )
	plot( as.mcmc.list(subPops(res,2), minPopLength=10), smooth=FALSE )
	resB <- thin(res,start=20)	
	#plot(res[[2]]$temp)
	
	#iPop=2
	.nGenThinned=rep(.nGen,.nPop)%/%.thin+1
	for( iPop in seq_along(res$pops) ){
		resPop <- res$pops[[iPop]]
		checkEquals( .nGenThinned[iPop], nrow(resPop$parms), msg="number of states in parms mismatch" )
		checkEquals( .nGenThinned[iPop], length(resPop$temp), msg="number of states in temp mismatch" )
		checkEquals( .nGenThinned[iPop], nrow(resPop$pAccept), msg="number of states in pAccept mismatch" )
		checkEquals( .nGenThinned[iPop], nrow(resPop$resLogDen), msg="number of states in resLogDen mismatch" )
		checkEquals( .nGenThinned[iPop], nrow(resPop$logDen), msg="number of states in logDen mismatch" )
		ZinitI <- pops0[[iPop]]$parms
		checkEquals( ZinitI[nrow(ZinitI),,], resPop$parms[1,,], msg="first row of parms should correspond to the last row of parms" )		
		checkEquals( 1, resPop$temp[.nGenThinned[iPop] ], , msg="end temperature need to be 1" )		
	}
	checkEquals( 1, res$pops[[1]]$temp[1], "first row of temperature of first population must be 1" )
	checkEquals( pops0[[2]]$T0, res$pops[[2]]$temp[1], "first row of temperature of second population must correspond correspond to prescribed argument" )
	
	.tmp.f <- function(){
		pop <- resB$pops[[1]]
		pop <- resB$pops[[2]]  # with higher temp
		matplot( pop$logDen[,1,], type="l" )	# acceptance rates of first block
		matplot( pop$pAccept[,1,], type="l" )	# acceptance rates of first block
		matplot( pop$temp, type="l" )	# temperature
		matplot( pop$resLogDen[,1,], type="l" )	# logLik obs
		matplot( pop$resLogDen[,2,], type="l" )	# logLik parms
		#require(twMiscRgl)
		plot( pop$parms[,"a",], pop$parms[,"b",], col=rainbow(.nGenThinned)  )
		plot( pop$parms[,"a",], pop$parms[,"b",], col=rainbow(dim(pop$parms)[3])[rep(1:dim(pop$parms)[3], each=dim(pop$parms)[1])]  )
		plot( density(pop$parms[,"a",]) )		
		plot( density(pop$parms[,"b",]) )
	}
	
	
	#mtrace(sfRemoteWrapper)
	#mtrace(.doDEMCStep)
	#mtrace(logDenGaussian)
	#mtrace(.doDEMCSteps )
	#mtrace(twDEMCInt)
	res <-  twDEMC( Zinit, nGen=.nGen, 
		fLogDen=logDenGaussian, argsFLogDen=argsFLogDen,
		nPop=.nPop
		,controlTwDEMC=list(thin=8)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	)
	#str(res)
	checkEquals( 8, res$thin )
	checkEquals( (.nGen%/%res$thin)+1,nrow(res$rLogDen))

	#windows(record=TRUE)
	rescoda <- as.mcmc.list(resB) 
	plot(rescoda, smooth=FALSE)
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
	
	.tmp.f <- function(){
		# with blocked updates chains do not sample the distribution
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
}

profile.f <- function(){
	#same setup as test.distinctLogDen
	.nGen=100
	#.thin=4
	#.nGen=3
	#mtrace(.checkBlock)
	#mtrace(.updateBlockTwDEMC)
	#mtrace(.updateBlocksTwDEMC)
	#mtrace(twDEMCBlockInt)
	#mtrace(.updateIntervalTwDEMCPar)
	res <- resAll <- twDEMCBlockInt( pops=pops, dInfos=dInfos, blocks=blocks, nGen=.nGen, controlTwDEMC=list(thin=.thin) )
	prof <- profr(
		res <- resAll <- twDEMCBlockInt( pops=pops, dInfos=dInfos, blocks=blocks, nGen=1, controlTwDEMC=list(thin=1), debugSequential=TRUE )
		,interval=0.001	)
	#plot(prof)
	prof2 <- subset(prof, start >= 0.01)
	plot(prof2, minlabel=0.05, angle=10)
	prof <- profr(
		res <- resAll <- twDEMCBlockInt( pops=pops, dInfos=dInfos, blocks=blocks, nGen=1, controlTwDEMC=list(thin=1), debugSequential=FALSE )
		,interval=0.001	)
	plot(prof)
	
	Rprof(interval=0.001)
	res <- resAll <- twDEMCBlockInt( pops=pops, dInfos=dInfos, blocks=blocks, nGen=.nGen, controlTwDEMC=list(thin=1)
	, debugSequential=TRUE )	
	Rprof(NULL)
	head(summaryRprof()$by.self)
	
	Rprof(interval=0.001)
	res <- resAll <- twDEMCBlockInt( pops=pops, dInfos=dInfos, blocks=blocks, nGen=.nGen, controlTwDEMC=list(thin=1)
	, debugSequential=FALSE )	
	Rprof(NULL)
	head(summaryRprof()$by.self)
	
	tmp <- system.time( res <- resAll <- twDEMCBlockInt( pops=pops, dInfos=dInfos, blocks=blocks, nGen=.nGen, controlTwDEMC=list(thin=4)) )
	# about 4 seconds for 200 runs 
	tmp2 <- system.time( res <- resAll <- twDEMCBlockInt( pops=pops, dInfos=dInfos, blocks=blocks, nGen=.nGen, controlTwDEMC=list(thin=4), debugSequential=FALSE) )
	# about 4 seconds for 200 runs with 4 cpus
	# overhead for one update is so about 20ms
}

test.extendRun <- function(){
	lmDummy <- lm( obs ~ xval, weights=1/sdObs^2)		# results without priors
	(.expTheta <- coef(lmDummy))
	(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
	(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=c("a","b")) )
	
	# nice starting values
	.nPop=4
	argsFLogDen <- list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,			### vector of data to compare with
		invCovar=invCovar,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
		#thetaPrior = thetaTrue,	### the prior estimate of the parameters
		#invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
		xval=xval
	)
	do.call( logDenGaussian, c(list(theta=theta0),argsFLogDen))
	
	Zinit <- initZtwDEMCNormal( theta0, .expCovTheta, nChainPop=4, nPop=.nPop)
	#dim(Zinit)
	
	.nGen=c(16,64,0,0)
	#mtrace(logDenGaussian)
	#mtrace(twDEMCBlockInt)
	res <-  concatPops( resBlock <- twDEMCBlock( Zinit, nGen=.nGen, 
			dInfos=list(
				d1 = list(fLogDen=logDenGaussian, argsFLogDen=argsFLogDen)
				,d2 = list(fLogDen=logDenGaussian, argsFLogDen=argsFLogDen)
			)
			,blocks = list(
				b1 = list(dInfoPos="d1", compPos="a")
				,b1 = list(dInfoPos="d1", compPos="b")
			)
			,nPop=.nPop
			,controlTwDEMC=list(thin=4)		
			,debugSequential=TRUE
			,doRecordProposals=TRUE
		), minPopLength=2)	# omitting the zero generation chain from  
	#str(res)
	rescoda <- as.mcmc.list(res) 
	getNSamples(resBlock)
	getNSamples(res)
	plot(rescoda, smooth=FALSE)
	#gelman.diag(rescoda)
	#summary(rescoda)

	checkEquals(.nGen, getNGen(resBlock) )
	checkEquals(.nGen, sapply(resBlock$pops, function(pop){ nrow(pop$Y)}) )
	
	.nGen2 <- c(0,64+32,16,0)	# will produce warnings for those populations 3,4 with 0 generations
	#mtrace(twDEMCBlockInt)
	resBlock2 <- twDEMCBlock(resBlock, nGen=.nGen2
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	)
	getNSamples(resBlock2)
	checkEquals( .nGen+.nGen2, getNGen(resBlock2) )
	checkEquals( .nGen+.nGen2, sapply( resBlock2$pops, function(pop){ nrow(pop$Y)}) )
	
	#mtrace(twDEMCBlock.twDEMCPops)
	resBlock2b <- twDEMCBlock(resBlock, nGen=.nGen2
		,debugSequential=TRUE
		#,doRecordProposals=TRUE
	)
	checkEquals( .nGen+.nGen2, getNGen(resBlock2b) )
	checkEquals( pmin(.nGen+.nGen2,128), sapply( resBlock2b$pops, function(pop){ nrow(pop$Y)}) )
}





