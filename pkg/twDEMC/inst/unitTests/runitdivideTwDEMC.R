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


test.runDivideTwDEMC <- function(){
	res <- divideTwDEMC(den2dCorDivideTwDEMC$sample, nGen=512, fLogDen=den2dCor )
	plot( b ~ a, as.data.frame(res[[1]]$sample), xlim=c(-0.5,2), ylim=c(-20,40) )
	points(0.8, 0, col="red" )	# theoretical maximum
	ss <- lapply(res , "[[" ,"sample" )	
	checkEquals( length(ss), 2, msg="wrong number of populations")	# two populations
	#ssPop <- ss[[1]]
	#ssPop <- ss[[2]]
	sapply(ss, function(ssPop){ 
		checkEquals( nrow(ssPop), 500/4, msg=" wrong number of generations in subPopulation") # number of generations 
		checkEquals( ncol(ssPop), 3, msg=" wrong number of parameters in subPopulation") 	 # number of parameters + rLogLik
		imax <- which.max(ssPop[,1])
		thetaHat <- ssPop[imax,]			
		checkMagnitude(thetaHat["a"], 0.8 ,msg=" wrong thetaHat[a]")
		checkInterval(thetaHat["b"], -5, +5,msg=" wrong thetaHat[b]")
		means <- apply( ssPop,2, mean)
		checkInterval( means["a"], 0, 1.2,msg=" wrong mean(a)")
	})	
}

test.initialShift <- function(){
	# start on a worse sample
	ss0 <- stackChains(den2dCorTwDEMC)
	#mtrace(divideTwDEMC)
	res <- divideTwDEMC(stackChainsPop(den2dCorTwDEMC), nGen=256, fLogDen=den2dCor )
	ssImpPops1 <- ssImpPops <- abind( lapply( res, "[[", "sample"), rev.along=0 )
	plot(density( ssImpPops[,"a",1]));lines(density( ssImpPops[,"a",2]),col="green"); lines(density( ss0[,"a"]),col="blue")
	checkTrue( all(apply(ssImpPops1[,"a",], 2, mean) > -0.8), msg="did not shift means(a) towards narrow part." )
}

