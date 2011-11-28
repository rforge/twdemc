.setUp <- function(){
	data(twLinreg1)
	data(twdemcEx1)
	attach(twLinreg1)
}

.tearDown <- function(){
	detach()
}

test.SubChainsAndcombinePops <- function(){
	#make sure twdemcEx1 has 2 populations of 4 chains
	ex1c <- concatPops(twdemcEx1)
	
	.nPops<-2
	#all subChains
	.pops <- lapply( 1:.nPops, function(iPop){ subChains(ex1c, iPop=iPop) })
	
	#same as specifying all chains
	#mtrace(subChains.twDEMC)
	.tmp <- subChains(ex1c,1:4)
	#lapply( names(.tmp), function(i){ print(i); checkEquals(.tmp[[i]], .pops[[1]][[i]]) })
	checkEquals(.pops[[1]],.tmp)
	
	#combine it again
	#mtrace(combinePops.twDEMC)
	#combinePops.twDEMC( .pops[[1]], .pops[[2]] )
	.combPop <- do.call( combinePops, .pops)
	.tmpexp<-ex1c;attr(.tmpexp,"batchCall")<-NULL	#batchCall only in original
	#for( i in  names(.tmpexp) ){print(i); checkEquals(.tmpexp[[i]], .combPop[[i]])	}
	#i="parms"
	#str(.tmpexp[[i]])
	#str(.combPop[[i]])
	checkEquals(.tmpexp, .combPop)
	
	#thinning
	#depends on subset
	.newThin=2*ex1c$thin
	.popThinned <- thin(.pops[[1]], newThin=.newThin)
	checkEquals(.newThin, .popThinned$thin)
	checkEquals(floor(100/.newThin)+1, nrow(.popThinned$logDen) )
	
	#checking the exception for different thinning attributes
	checkException( do.call( combinePops, c(.pops, list(.popThinned)) ))
}


test.stackChains <- function(){
	ex1c <- concatPops(twdemcEx1)
	.tmp <- stackChains(ex1c)
	checkTrue( is.matrix(.tmp) )
	checkEquals( colnames(ex1c$parms), colnames(.tmp[,-attr(.tmp,"nBlock")]) )
	checkEquals( prod( dim(ex1c$logDen) ), nrow(.tmp) )
}

test.thin.twDEMC <- function(){
	# test thin.twDEMC
	#mtrace(thin.twDEMC)
	# test modifying starting value (cutting burnin)
	ex1c <- concatPops(twdemcEx1)
	ex1c$nGenBurnin=8*ex1c$thin
	thinned <- thin(ex1c, start=ex1c$nGenBurnin)	# removing burnin period
	checkEquals(c(getNGen(ex1c)-ex1c$nGenBurnin,ex1c$thin,getNSamples(ex1c)-8,0)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	# test giving time a bit between thinning intervals
	thinned <- thin(ex1c, start=1)	# remove first thinning interval (because not starting from time zero)
	checkEquals(c(getNGen(ex1c)-ex1c$thin,ex1c$thin,getNSamples(ex1c)-1, ex1c$nGenBurnin-ex1c$thin)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	# test modifying end value (cutting from the end of the chain)
	thinned <- thin(ex1c, end=getNGen(ex1c))	# no change
	checkEquals(c(getNGen(ex1c)-0,ex1c$thin,getNSamples(ex1c)-0,ex1c$nGenBurnin-0)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	thinned <- thin(ex1c, end=getNGen(ex1c)+5)	# end behind sample: no change
	checkEquals(c(getNGen(ex1c)-0,ex1c$thin,getNSamples(ex1c)-0,ex1c$nGenBurnin-0)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	thinned <- thin(ex1c, end=getNGen(ex1c)-2)	# remove entire last thinning interval
	checkEquals(c(getNGen(ex1c)-1*ex1c$thin,ex1c$thin,getNSamples(ex1c)-1,ex1c$nGenBurnin-0)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	thinned <- tM2 <- thin(ex1c, start=2, end=getNGen(ex1c)-2)	# starting from 2: removing first and last thinning interval
	checkEquals(c(getNGen(ex1c)-2*ex1c$thin,ex1c$thin,getNSamples(ex1c)-2,ex1c$nGenBurnin-1*ex1c$thin)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	# test new thinning interval
	.newThin=2*ex1c$thin
	thinned <- t8 <- thin(ex1c, newThin=.newThin)	
	checkEquals(c( floor(getNGen(ex1c)/.newThin)*.newThin,.newThin,(getNSamples(ex1c)-1)%/%2+1,ex1c$nGenBurnin)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	thinned <- thin(ex1c, newThin=.newThin, start=2, end=getNGen(ex1c)-2)	
	checkEquals(c( floor(getNGen(tM2)/.newThin)*.newThin,.newThin,(getNSamples(tM2)-1)%/%2+1,tM2$nGenBurnin)
		, c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin[1] )
		, checkNames=FALSE )
	
}

#XXTODO
.tmp.f <- function(){
popMeansTwDEMC
}



