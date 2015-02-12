.tmp.f <- function(){
    twUtestF("twCalcLogDenPar") 
}

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
    ex1c <- concatPops(twdemcEx1)
    xProp <- stackChains(ex1c$parms)
    #trace(twCalcLogDenPar, recover)    #untrace(twCalcLogDenPar)
    # without SSResid, without intermediate
    # two costs, several rows
    str(res <- twCalcLogDenPar(function(x){2*x},xProp,debugSequential=TRUE))
    checkEquals( dim(res$logDenComp), dim(xProp) )
    checkEquals( as.vector(res$logDenComp), as.vector(xProp*2) )
    # two costs, one row
    str(res <- twCalcLogDenPar(function(x){2*x},xProp[1, ,drop=FALSE],debugSequential=TRUE))
    checkEquals( dim(res$logDenComp), c(1,2) )
    checkEquals( as.vector(res$logDenComp), as.vector(xProp[1,]*2) )
    # one cost, several rows
    str(res <- twCalcLogDenPar(function(x){2*x[1]},xProp,debugSequential=TRUE))
    checkEquals( dim(res$logDenComp), c(nrow(xProp),1) )
    checkEquals( as.vector(res$logDenComp), as.vector(xProp[,1]*2) )
    # one cost, one row
    str(res <- twCalcLogDenPar(function(x){2*x[1]},xProp[1,,drop=FALSE],debugSequential=TRUE))
    checkEquals( dim(res$logDenComp), c(1,1) )
    checkEquals( as.vector(res$logDenComp), as.vector(xProp[1,1]*2) )

    # with logDenCompX, vector result, providing only truncated set of parameters
    .logDenCompX0 <- matrix(1:nrow(xProp), nrow=nrow(xProp), dimnames=list(steps=NULL,"a"))
    .logDenCompX0[1:5, ,drop=FALSE]
    #function must take two arguments
    #mtrace(twCalcLogDenPar)
    str(res <- twCalcLogDenPar(function(y1,logDenCompAcc){y1["a"]=y1["a"]*logDenCompAcc["a"]; y1},xProp,logDenCompX=.logDenCompX0, intResComp=colnames(.logDenCompX0), debugSequential=TRUE))
    checkEquals( nrow(xProp), length(res$logDen) )
    .exp <- xProp; .exp[,"a"] = xProp[,"a"]*.logDenCompX0[,"a"]
    names(dimnames(.exp)) <- NULL
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

test.twCalcLogDensPars_namedComp <- function(){
    # testing named or unnmaed return components
    lDen1NoNames <- function(theta)( 1 )
    lDen2NoNames <- function(theta)( c(1,2) )
    lDen1Names <- function(theta)( c(r1=1) )
    lDen2Names <- function(theta)( c(r1=1,r2=2) )
    
    xProp1 <- matrix(1, nrow=1, ncol=1)
    xProp2 <- matrix(1, nrow=2, ncol=1)

    res <- twCalcLogDensPar(dInfos=list(list(fLogDen=lDen1NoNames)), xProp1)
    checkEquals( "d1_1", colnames(res$logDenComp) )
    res <- twCalcLogDensPar(dInfos=list(list(fLogDen=lDen1NoNames)), xProp2)
    checkEquals( "d1_1", colnames(res$logDenComp) )
    res <- twCalcLogDensPar(dInfos=list(list(fLogDen=lDen2NoNames)), xProp1)
    checkEquals( c("d1_1","d1_2"), colnames(res$logDenComp) )
    res <- twCalcLogDensPar(dInfos=list(list(fLogDen=lDen2NoNames)), xProp2)
    checkEquals( c("d1_1","d1_2"), colnames(res$logDenComp) )
    
    res <- twCalcLogDensPar(dInfos=list(list(fLogDen=lDen1Names)), xProp1)
    checkEquals( "r1", colnames(res$logDenComp) )
    res <- twCalcLogDensPar(dInfos=list(list(fLogDen=lDen1Names)), xProp2)
    checkEquals( "r1", colnames(res$logDenComp) )
    res <- twCalcLogDensPar(dInfos=list(list(fLogDen=lDen2Names)), xProp1)
    checkEquals( c("r1","r2"), colnames(res$logDenComp) )
    res <- twCalcLogDensPar(dInfos=list(list(fLogDen=lDen2Names)), xProp2)
    checkEquals( c("r1","r2"), colnames(res$logDenComp) )
    
    res <- twCalcLogDensPar(dInfos=list(den1=list(fLogDen=lDen1NoNames)), xProp1)
    checkEquals( "den1_1", colnames(res$logDenComp) )
    res <- twCalcLogDensPar(dInfos=list(den1=list(fLogDen=lDen1Names)), xProp1)
    checkEquals( "r1", colnames(res$logDenComp) )

    res <- twCalcLogDensPar(dInfos=list(
                 den1=list(fLogDen=lDen2Names)
                ,den2=list(fLogDen=lDen2NoNames)
                ,den3=list(fLogDen=lDen1NoNames)
    ), xProp1)
    checkEquals( c("r1","r2","den2_1","den2_2","den3_1"), colnames(res$logDenComp) )
    
}


