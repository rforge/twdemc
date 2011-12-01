.setUp <- function(){
	data(twLinreg1)
	attach( twLinreg1 )
}

.tearDown <- function(){
	detach( twLinreg1 )
}

test.initZtwDEMCNormal <- function(){
	.nChainPop=4
	.nPop=2
	.nPar=length(theta0)
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChainPop=.nChainPop, nPop=.nPop)
	head(Zinit[,,1])
	checkEquals( c(calcM0twDEMC(.nPar,.nChainPop), .nPar, .nChainPop*.nPop), dim(Zinit) )
}

test.initZtwDEMCNormalPar1 <- function(){
	# test for parameter vector of length 1
	.nChainPop=4
	.nPop=2
	.nPar=1
	Zinit <- initZtwDEMCNormal( theta0[1], diag(sdTheta^2)[1], nChainPop=.nChainPop, nPop=.nPop)
	Zinit[,,1,drop=FALSE]
	checkEquals( c(calcM0twDEMC(.nPar,.nChainPop), .nPar, .nChainPop*.nPop), dim(Zinit) )
}

test.initZtwDEMCNormalPar1Chain1 <- function(){
	# test for parameter vector of length 1 and only one chain
	.nChainPop=1
	.nPop=1
	.nPar=1
	Zinit <- initZtwDEMCNormal( theta0[1], diag(sdTheta^2)[1], nChainPop=.nChainPop, nPop=.nPop)
	Zinit[,,1,drop=FALSE]
	checkEquals( c(calcM0twDEMC(.nPar,.nChainPop), .nPar, .nChainPop*.nPop), dim(Zinit) )
}


#mtrace(initZtwDEMCSub.array)
test.initZtwDEMCSub <- function(){
	data(twdemcEx1)
	#x <- concatPops(twdemcEx1)
	#mtrace(stackChains.twDEMC)
	#s1 <- stackChains(x)
	#mtrace(initZtwDEMCSub.twDEMC)
	#mtrace(initZtwDEMCSub.matrix)
	Zinit <- initZtwDEMCSub( concatPops(twdemcEx1), vars=c("a","b") )
	checkEquals( c(4,2,8), dim(Zinit) )
}

test.initZExt <- function(){
	data(twdemcEx1)
	#mtrace(initZtwDEMCExt.matrix)
	Zinita <- initZtwDEMCSub( concatPops(twdemcEx1), vars=c("a") )
	checkEquals( c(4,1,8), dim(Zinita) )
	
	#now extend the subsample again by subsampling normal b
	Zinit <- initZtwDEMCExt( stackChains(Zinita), thetaPrior=theta0, covarTheta=diag(sdTheta^2),nChains=dim(Zinita)[3], nPop=getNPops(twdemcEx1)) 
	checkEquals( c(4,2,8), dim(Zinit) )
}

test.constrainNStack <- function(){
	pss <- stackChains(concatPops(twdemcEx1))
	#normpoptBest <-  twExtractFromLastDims(twdemcEx1$parms, which.max( twdemcEx1$logDen) )[,1]
	normpoptBest <-  pss[ which.max( pss[,1]), -1]
	pss2 <- constrainNStack(pss, thetaPrior=normpoptBest[c("a","b")], n = 80 )
	checkEquals( c(80,3), dim(pss2) )
	
	pss2 <- constrainNStack(pss, thetaPrior=normpoptBest[c("a")], n = 80 )
	checkEquals( c(80,3), dim(pss2) )

	#mtrace(constrainNStack)
	pss2 <- constrainNStack(pss, thetaPrior=c(), n = 80, returnAlpha=TRUE ) # just subsampling
	checkEquals( c(80,3), dim(pss2$res) )
}

test.constrainCfStack <- function(){
	pss <- stackChains(concatPops(twdemcEx1))
	normpoptBest <-  pss[ which.max( pss[,1]), -1]
	
	.alpha=0.95
	expN <- .alpha*nrow(pss)
	
	#mtrace(constrainCfStack)
	pss2 <- constrainCfStack(pss, thetaPrior=normpoptBest[c("a","b")], alpha=.alpha )
	checkMagnitude( expN, nrow(pss2) )
	checkTrue( nrow(pss2) <= ceiling(expN))
	
	pss2 <- constrainCfStack(pss, thetaPrior=normpoptBest[c("a")], alpha=.alpha)
	checkMagnitude( expN, nrow(pss2) )
	checkTrue( nrow(pss2) <= ceiling(expN))
	plot(density(pss[,"a"]))
	lines(density(pss2[,"a"]),col="red")
	
	#mtrace(constrainNStack)
	pss2 <- constrainCfStack(pss, thetaPrior=c(), alpha=.alpha )
	checkEquals( pss, pss2 )
}