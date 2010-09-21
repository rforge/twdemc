#twUtest(replaceZinitNonFiniteLogLiks)

.setUp <- function(){
	data(twdemcEx1)
}

.tearDown <- function(){
}

test.default <- function(){
	.Zinit <- twdemcEx1$parms
	.rLogLik <- .rLogLikNA <- twdemcEx1$rLogLik
	.naCases <- 1:10
	.naChains <- c(1,8)
	.rLogLikNA[.naCases,.naChains] <- NA 
	.Zinit2 <- replaceZinitNonFiniteLogLiks(.Zinit,.rLogLikNA)
	checkEquals( dim(.Zinit), dim(.Zinit2))
	checkTrue( all(.Zinit2[,.naCases,.naChains] != .Zinit[,.naCases,.naChains]) ) 
	checkTrue( all(.Zinit2[,-.naCases,.naChains] == .Zinit[,-.naCases,.naChains]) ) 
	checkTrue( all(.Zinit2[,,-.naChains] == .Zinit[,,-.naChains]) ) 
}

test.lastRow <- function(){
	.Zinit <- array( 1:20, dim=c(1,5,4) )
	.Zinit[,5,1:2] <- NA	#first two chains of last row
	.Zinit[,2:3,] <- Inf	#second and third state of all chains, leaving 1 and 4th state finite		
	fLogLik <- function(x){ x }
	checkTrue( !is.finite(fLogLik( .Zinit[,5,1])) )
	checkTrue( is.finite(fLogLik( .Zinit[,5,3])) )
	checkTrue( !is.finite(fLogLik( .Zinit[,2,3])) )
	
	#mtrace(replaceZinitNonFiniteLogLiksLastStep)
	res<-NULL; res <- replaceZinitNonFiniteLogLiksLastStep(.Zinit,fLogLik,nPops=2)
	checkEquals( dim(res$Zinit), dim(res$Zinit))
	rLogLik <- twCalcLogLikPar( fLogLik, t(adrop(res$Zinit[,5,,drop=FALSE],2)) )$logLik
	checkTrue( all(is.finite(rLogLik)))
	checkEquals( length(rLogLik), length(unique(rLogLik)) )	#distinct states
	checkEquals( rLogLik, res$rLogLik )
}

.tmp.f <- function(){
	mtrace(replaceZinitNonFiniteLogLiks)
	mtrace(replaceZinitCases)
}

