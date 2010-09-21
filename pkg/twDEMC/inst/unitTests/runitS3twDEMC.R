.setUp <- function(){
	data(twLinreg1)
	attach(twLinreg1)
}

.tearDown <- function(){
	detach()
}



test.combinePops <- function(){
	#make sure twdemcEx1 has 2 populations of 4 chains
	
	.nPops<-2
	#all subChains
	.pops <- lapply( 1:.nPops, function(iPop){ subChains(twdemcEx1, iPop=iPop) })
	
	#same as specifying all chains
	#mtrace(subChains.twDEMC)
	.tmp <- subChains(twdemcEx1,1:4)
	checkEquals(.pops[[1]],.tmp)
	
	#combine it again
	#mtrace(combinePops.twDEMC)
	.combPop <- do.call( combinePops, .pops)
	.tmpexp<-twdemcEx1;attr(.tmpexp,"batchCall")<-NULL	#batchCall only in original
	checkEquals(.tmpexp, .combPop)
	
	#thinning
	.newThin=10
	.popThinned <- thin(.pops[[1]], newThin=.newThin)
	checkEquals(.newThin, .popThinned$thin)
	checkEquals(100/.newThin+1, nrow(.popThinned$rLogLik) )
	
	#checking the exception for different thinning attributes
	checkException( do.call( combinePops, c(.pops, list(.popThinned)) ))
}

test.stackChains <- function(){
	.tmp <- stackChains(twdemcEx1)
	checkTrue( is.matrix(.tmp) )
	checkEquals( c("rLogLik",rownames(twdemcEx1$parms)), colnames(.tmp) )
	checkEquals( prod( dim(twdemcEx1$rLogLik) ), nrow(.tmp) )
}



