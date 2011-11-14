.setUp <- function(){
	data(twLinreg1)
	data(twdemcEx1)
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
	#combinePops.twDEMC( .pops[[1]], .pops[[2]] )
	.combPop <- do.call( combinePops, .pops)
	.tmpexp<-twdemcEx1;attr(.tmpexp,"batchCall")<-NULL	#batchCall only in original
	checkEquals(.tmpexp, .combPop)
	
	#thinning
	.newThin=10
	.popThinned <- thin(.pops[[1]], newThin=.newThin)
	checkEquals(.newThin, .popThinned$thin)
	checkEquals(100/.newThin+1, nrow(.popThinned$rLogDen) )
	
	#checking the exception for different thinning attributes
	checkException( do.call( combinePops, c(.pops, list(.popThinned)) ))
}

test.stackChains <- function(){
	.tmp <- stackChains(twdemcEx1)
	checkTrue( is.matrix(.tmp) )
	checkEquals( c("rLogDen",rownames(twdemcEx1$parms)), colnames(.tmp) )
	checkEquals( prod( dim(twdemcEx1$rLogDen) ), nrow(.tmp) )
}

test.thin.twDEMC <- function(){
	# test thin.twDEMC
	data(twdemcEx1)
	#mtrace(thin.twDEMC)
	# test modifying starting value (cutting burnin)
	thinned <- thin(twdemcEx1, start=twdemcEx1$nGenBurnin)	# removing burnin period
	checkEquals(c(70,5,15,0)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	# test giving time a bit between thinning intervals
	thinned <- thin(twdemcEx1, start=1)	# remove first thinning interval (because not starting from time zero)
	checkEquals(c(95,5,20,25)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	# test modifying end value (cutting from the end of the chain)
	thinned <- thin(twdemcEx1, end=getNGen(twdemcEx1))	# no change
	checkEquals(c(100,5,21,30)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	thinned <- thin(twdemcEx1, end=getNGen(twdemcEx1)+5)	# end behind sample: no change
	checkEquals(c(100,5,21,30)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	thinned <- thin(twdemcEx1, end=getNGen(twdemcEx1)-2)	# remove entire last thinning interval
	checkEquals(c(95,5,20,30)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	thinned <- thin(twdemcEx1, start=2, end=getNGen(twdemcEx1)-2)	# remove entire last thinning interval
	checkEquals(c(90,5,19,25)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	# test new thinning interval
	thinned <- thin(twdemcEx1, newThin=10)	
	checkEquals(c(100,10,11,30)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	thinned <- thin(twdemcEx1, newThin=10, start=2, end=getNGen(twdemcEx1)-2)	
	checkEquals(c(80,10,9,20)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	
}



