#twUtest(replaceZinitNonFiniteLogDens)

.setUp <- function(){
	data(twdemcEx1)
}

.tearDown <- function(){
}

test.default <- function(){
	.Zinit <- twdemcEx1$parms
	.rLogDen <- .rLogDenNA <- twdemcEx1$rLogDen
	.naCases <- 1:10
	.naChains <- c(1,8)
	.rLogDenNA[.naCases,.naChains] <- NA 
	.Zinit2 <- replaceZinitNonFiniteLogDens(.Zinit,.rLogDenNA)
	checkEquals( dim(.Zinit), dim(.Zinit2))
	checkTrue( all(.Zinit2[,.naCases,.naChains] != .Zinit[,.naCases,.naChains]) ) 
	checkTrue( all(.Zinit2[,-.naCases,.naChains] == .Zinit[,-.naCases,.naChains]) ) 
	checkTrue( all(.Zinit2[,,-.naChains] == .Zinit[,,-.naChains]) ) 
}

test.lastRow <- function(){
	.Zinit <- array( 1:20, dim=c(1,5,4) )
	.Zinit[,5,1:2] <- NA	#first two chains of last row
	.Zinit[,2:3,] <- Inf	#second and third state of all chains, leaving 1 and 4th state finite		
	fLogDen <- function(x){ x }
	checkTrue( !is.finite(fLogDen( .Zinit[,5,1])) )
	checkTrue( is.finite(fLogDen( .Zinit[,5,3])) )
	checkTrue( !is.finite(fLogDen( .Zinit[,2,3])) )
	
	#mtrace(replaceZinitNonFiniteLogDensLastStep)
	res<-NULL; res <- replaceZinitNonFiniteLogDensLastStep(.Zinit,fLogDen,nPops=2)
	checkEquals( dim(res$Zinit), dim(res$Zinit))
	rLogDen <- twCalcLogDenPar( fLogDen, t(adrop(res$Zinit[,5,,drop=FALSE],2)) )$logDen
	checkTrue( all(is.finite(rLogDen)))
	checkEquals( length(rLogDen), length(unique(rLogDen)) )	#distinct states
	checkEquals( rLogDen, res$rLogDen )
}

.tmp.f <- function(){
	mtrace(replaceZinitNonFiniteLogDens)
	mtrace(replaceZinitCases)
}

