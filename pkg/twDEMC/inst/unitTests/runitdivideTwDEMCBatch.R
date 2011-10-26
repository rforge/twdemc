.setUp <- function(){
	data(den2dCorDivideTwDEMC)
	data(den2dCorTwDEMC)
	#attach(twLinreg1)
	#attach(twdemcEx1)
}

.tearDown <- function(){
	#detach( twdemcEx1 )
	#detach( twLinreg1 )
}


test.runDivideTwDEMCBatch <- function(){
	res <- divideTwDEMCBatch(den2dCorDivideTwDEMC$sample, nGen=512*2, fLogDen=den2dCor )
	plot( b ~ a, as.data.frame(res$sample[,,1]), xlim=c(-0.5,2), ylim=c(-20,40) )
	points(0.8, 0, col="red" )	# theoretical maximum
	#x <- A <- res$sample
	ss <- twListArrDim(res$sample)	
	checkEquals( length(ss), 2, msg="wrong number of populations")	# two populations
	#ssPop <- ss[[1]]
	sapply(ss, function(ssPop){ 
		checkTrue( nrow(ssPop) > 500/4, msg="wrong number of generations in subPopulation") # number of generations 
		checkEquals( ncol(ssPop), 3, msg="wrong number of parameters in subPopulation") 	 # number of parameters + rLogLik
		imax <- which.max(ssPop[,1])
		thetaHat <- ssPop[imax,]			
		checkMagnitude(thetaHat["a"], 0.8 ,msg="wrong thetaHat[a]")
		checkInterval(thetaHat["b"], -5, +5,msg="wrong thetaHat[b]")
		means <- apply( ssPop,2, mean)
		checkInterval( means["a"], 0, 1.2,msg="wrong mean(a)")
	})	
}
