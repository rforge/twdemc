.setUp <- function(){
	#data(twLinreg1)
	#attach(twLinreg1)
	#attach(twdemcEx1)
}

.tearDown <- function(){
	#detach( twdemcEx1 )
	#detach( twLinreg1 )
}

test.twCalcLogDenPar <- function(){
    data(twdemcEx1)
    xProp <- stackChains(concatPops(twdemcEx1)$parms)
	#mtrace(twCalcLogDenPar)
	#logDenCompXStacked <- stackChains(logDenCompX)
	str(res <- twCalcLogDenPar(function(x){2*x},xProp,debugSequential=TRUE))
	
	# with logDenCompX, vector result, providing only truncated set of parameters
	.logDenCompX0 <- matrix(1:nrow(xProp), nrow=nrow(xProp), dimnames=list(steps=NULL,"a"))
	.logDenCompX0[1:5,,drop=FALSE]
	#function must take two arguments, if provided previous results
	#mtrace(twCalcLogDenPar)
	str(res <- twCalcLogDenPar(
			function(y1,logDenCompAcc){	y1["a"]=y1["a"]*logDenCompAcc["a"]; y1}
			,xProp
			, logDenCompX=.logDenCompX0
			, intResComp=colnames(.logDenCompX0)
			, debugSequential=TRUE
	))

	checkEquals( nrow(xProp), length(res$logDen) )
	.exp <- xProp; .exp[,"a"] = xProp[,"a"]*.logDenCompX0[,"a"]; names(dimnames(.exp)) <- NULL
	checkEquals(.exp, res$logDenComp)
}

test.twCalcLogDensPars <- function(){
    data(twTwoDenEx1)
    thresholdCovar = 0.3	# the true value
    #thresholdCovar = 0		# the effective model that glosses over this threshold
    thetaMean <- twTwoDenEx1$thetaTrue
    sdTheta <- (thetaMean*0.3)
    #
    argsFLogDen <- list(
            thresholdCovar=thresholdCovar
            , thetaPrior=thetaMean
            , invCovarTheta=sdTheta^2
            , twTwoDenEx=twTwoDenEx1
    )            
    #
    .nPop=2
    .nChainPop=4
    ZinitPops <- with(twTwoDenEx1, initZtwDEMCNormal( thetaMean*1.5, diag(sdTheta^2), nChainPop=.nChainPop, nPop=.nPop))
    xProp <- stackChains(ZinitPops)
    dInfos=list(
            dSparse=list(fLogDen=denSparsePrior, argsFLogDen=argsFLogDen, intermediate="modTwoDenEx")
            ,dRich=list(fLogDen=denRichPrior, argsFLogDen=argsFLogDen, intermediate="modTwoDenEx")
            ,dSparse2=list(fLogDen=denSparsePrior, argsFLogDen=argsFLogDen)   # test no intermediate specified
    )
    res <- twCalcLogDensPar(dInfos, xProp)
    checkEquals( nrow(xProp), nrow(res$logDenComp),  )
    checkEquals( 6, ncol(res$logDenComp) )
    
    res <- twCalcLogDensPar(dInfos, xProp[1, ,drop=FALSE])   # test on degenerate matrix
    checkEquals( 1, nrow(res$logDenComp) )
    checkEquals( 6, ncol(res$logDenComp) )
}

