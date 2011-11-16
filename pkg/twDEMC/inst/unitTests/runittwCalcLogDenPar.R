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

test.twCalcLogDenPar <- function(){
	xProp <- stackChains(twdemcEx1$parms)
	#mtrace(twCalcLogDenPar)
	#logDenCompXStacked <- stackChains(logDenCompX)
	str(res <- twCalcLogDenPar(function(x){2*x},xProp,debugSequential=TRUE))
	
	# with logDenCompX, vector result, providing only truncated set of parameters
	.logDenCompX0 <- matrix(1:nrow(xProp), nrow=nrow(xProp), dimnames=list(steps=NULL,"a"))
	.logDenCompX0[1:5,,drop=FALSE]
	#function must take two arguments
	#mtrace(twCalcLogDenPar)
	str(res <- twCalcLogDenPar(function(y1,logDenCompAcc){y1["a"]=y1["a"]*logDenCompAcc["a"]; y1},xProp,logDenCompX=.logDenCompX0, intResCompNames=colnames(.logDenCompX0), debugSequential=TRUE))
	checkEquals( nrow(xProp), length(res$logDen) )
	.exp <- xProp; .exp[,"a"] = xProp[,"a"]*.logDenCompX0[,"a"]
	checkEquals(.exp, res$logDenComp)
}

