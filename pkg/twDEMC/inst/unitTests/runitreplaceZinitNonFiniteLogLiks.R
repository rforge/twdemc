#twUtest(replaceZinitNonFiniteLogDens)

.setUp <- function(){
	data(twdemcEx1)
}

.tearDown <- function(){
}

test.default <- function(){
	mc1 <- concatPops(twdemcEx1)
	.Zinit <- mc1$parms
	.logDen <- .logDenNA <- mc1$logDen
	.naCases <- 1:10		# from 26 in mc1
	.naChains <- c(1,8)
	.logDenNA[.naCases,1,.naChains] <- NA 
	.Zinit2 <- replaceZinitNonFiniteLogDens(.Zinit,.logDenNA)
	checkEquals( dim(.Zinit), dim(.Zinit2))
	checkTrue( all(.Zinit2[.naCases,,.naChains] != .Zinit[.naCases,,.naChains]) ) 
	checkTrue( all(.Zinit2[-.naCases,,.naChains] == .Zinit[-.naCases,,.naChains]) ) 
	checkTrue( all(.Zinit2[,,-.naChains] == .Zinit[,,-.naChains]) ) 
}

test.lastRow <- function(){
	.Zinit <- array( 1:20, dim=c(5,1,4) )
	.Zinit[5,,1:2] <- NA	#first two chains of last row
	.Zinit[2:3,,] <- Inf	#second and third state of all chains, leaving 1 and 4th state finite		
	fLogDen <- function(x){ x }
	checkTrue( !is.finite(fLogDen( .Zinit[5,,1])) )
	checkTrue( is.finite(fLogDen( .Zinit[5,,3])) )
	checkTrue( !is.finite(fLogDen( .Zinit[2,,3])) )
	
	#mtrace(replaceZinitNonFiniteLogDensLastStep)
	res<-NULL; res <- replaceZinitNonFiniteLogDensLastStep(.Zinit,fLogDen,nPop=2)
	checkEquals( dim(res$Zinit), dim(res$Zinit))
	rLogDen <- twCalcLogDenPar( fLogDen, t(adrop(res$Zinit[5,, ,drop=FALSE],2)) )$logDen
	checkTrue( all(is.finite(rLogDen)))
	checkEquals( length(rLogDen), length(unique(rLogDen)) )	#distinct states
	checkEquals( rLogDen, res$rLogDen )
}

.tmp.f <- function(){
	mtrace(replaceZinitNonFiniteLogDens)
	mtrace(replaceZinitCases)
}

