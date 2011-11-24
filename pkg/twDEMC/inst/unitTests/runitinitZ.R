.setUp <- function(){
	data(twLinreg1)
	attach( twLinreg1 )
}

.tearDown <- function(){
	detach( twLinreg1 )
}

test.initZtwDEMCNormal <- function(){
	.nChainsPop=4
	.nPops=2
	.nPar=length(theta0)
	Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChainsPop=.nChainsPop, nPops=.nPops)
	head(Zinit[,,1])
	checkEquals( c(calcM0twDEMC(.nPar,.nChainsPop), .nPar, .nChainsPop*.nPops), dim(Zinit) )
}

test.initZtwDEMCNormalPar1 <- function(){
	# test for parameter vector of length 1
	.nChainsPop=4
	.nPops=2
	.nPar=1
	Zinit <- initZtwDEMCNormal( theta0[1], diag(sdTheta^2)[1], nChainsPop=.nChainsPop, nPops=.nPops)
	Zinit[,,1,drop=FALSE]
	checkEquals( c(calcM0twDEMC(.nPar,.nChainsPop), .nPar, .nChainsPop*.nPops), dim(Zinit) )
}

test.initZtwDEMCNormalPar1Chain1 <- function(){
	# test for parameter vector of length 1 and only one chain
	.nChainsPop=1
	.nPops=1
	.nPar=1
	Zinit <- initZtwDEMCNormal( theta0[1], diag(sdTheta^2)[1], nChainsPop=.nChainsPop, nPops=.nPops)
	Zinit[,,1,drop=FALSE]
	checkEquals( c(calcM0twDEMC(.nPar,.nChainsPop), .nPar, .nChainsPop*.nPops), dim(Zinit) )
}


test.initZtwDEMCSub <- function(){
	data(twdemcEx1)
	Zinit <- initZtwDEMCSub( concatPops(twdemcEx1)$parms, c("a","b") )
	checkEquals( c(2,4,8), dim(Zinit) )
}

test.initZExt <- function(){
	data(twdemcEx1)
	#mtrace(initZtwDEMCExt.matrix)
	Zinita <- initZtwDEMCSub( twdemcEx1, c("a") )
	checkEquals( c(1,4,8), dim(Zinita) )
	
	#now extend the subsample again by subsampling normal b
	Zinit <- initZtwDEMCExt( stackChains(Zinita), theta0, diag(sdTheta^2),nChains=dim(Zinita)[3]) 
	checkEquals( c(2,4,8), dim(Zinit) )
}

test.constrainNStack <- function(){
	pss <- stackChains(twdemcEx1)
	normpoptBest <- normpoptBest <- twExtractFromLastDims(twdemcEx1$parms, which.max( twdemcEx1$rLogDen) )[,1]
	pss2 <- constrainNStack(pss, thetaPrior=normpoptBest[c("a","b")], n = 80 )
	checkEquals( c(80,3), dim(pss2) )
	
	pss2 <- constrainNStack(pss, thetaPrior=normpoptBest[c("a")], n = 80 )
	checkEquals( c(80,3), dim(pss2) )

	#mtrace(constrainNStack)
	pss2 <- constrainNStack(pss, thetaPrior=c(), n = 80, returnAlpha=TRUE )
	checkEquals( c(80,3), dim(pss2$res) )
}

test.constrainCfStack <- function(){
	pss <- stackChains(twdemcEx1)
	normpoptBest <- normpoptBest <- twExtractFromLastDims(twdemcEx1$parms, which.max( twdemcEx1$rLogDen) )[,1]
	
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