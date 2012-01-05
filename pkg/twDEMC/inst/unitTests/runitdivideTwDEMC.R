.setUp <- function(){
	data(den2dCorEx)
	#attach(twLinreg1)
	#attach(twdemcEx1)
}

.tearDown <- function(){
	#detach( twdemcEx1 )
	#detach( twLinreg1 )
}


test.divideTwDEMCStep <- function(){
	.nGen <- 256
	res <- divideTwDEMCStep(den2dCorEx$mcSubspaces0, nGen=.nGen, dInfos=list(list(fLogDen=den2dCor)), debugSequential=TRUE )
	mc1 <- stackPopsInSpace(res$resTwDEMC) 
	.tmp.f <- function(){
		plot( b ~ a, as.data.frame(stackChains(mc1$pops[[1]]$parms)), xlim=c(-0.5,2), ylim=c(-20,40) )
		points(0.8, 0, col="red" )	# theoretical maximum
		getNSamples(res$resTwDEMC)
	}
	checkEquals( getNSpaces(den2dCorEx$mcSubspaces0), getNPops(mc1), msg="wrong number of spaces" )
	
	iPop=1
	sapply(1:getNPops(mc1), function(iPop){
		pop <- mc1$pops[[iPop]]
		checkTrue( nrow(pop$parms) >  256/4-10, msg=" wrong number of generations in subPopulation") # number of generations 
		ss <- stackChains(subPops(mc1, iPops=iPop ))
		imax <- which.max(ss[,1])
		thetaHat <- ss[imax,]			
		checkInterval(thetaHat["a"], -0.8, +0.8,msg="wrong thetaHat[a]")
		checkInterval(thetaHat["b"], -20, +20,msg=" wrong thetaHat[b]")
		checkTrue( mean(ss[,"a"]) > -0.8, msg="did not shift means(a) towards narrow part." )
	})	
}
#twUtestF(getSubSpaces)

test.divideTwDEMCSteps <- function(){
	.nGenBatch <- 256
	res <- divideTwDEMCSteps(den2dCorEx$mcSubspaces0, nGen=.nGenBatch*2, nGenBatch=.nGenBatch, dInfos=list(list(fLogDen=den2dCor)), debugSequential=TRUE )
	mc1 <- stackPopsInSpace(res$resTwDEMC) 
	.tmp.f <- function(){
		plot( b ~ a, as.data.frame(stackChains(mc1$pops[[1]]$parms)), xlim=c(-0.5,2), ylim=c(-20,40) )
		points(0.8, 0, col="red" )	# theoretical maximum
		getNSamples(res$resTwDEMC)
	}
	checkEquals( getNSpaces(den2dCorEx$mcSubspaces0), getNPops(mc1), msg="wrong number of spaces" )
	
	iPop=1
	sapply(1:getNPops(mc1), function(iPop){
			#pop <- mc1$pops[[iPop]]
			#checkTrue( nrow(pop$parms) >  256/4-10, msg=" wrong number of generations in subPopulation") # number of generations 
			ss <- stackChains(subPops(mc1, iPops=iPop ))
			imax <- which.max(ss[,1])
			thetaHat <- ss[imax,]			
			checkMagnitude(thetaHat["a"], 0.8 ,msg=" wrong thetaHat[a]")
			checkInterval(thetaHat["b"], -5, +5,msg=" wrong thetaHat[b]")
			checkTrue( mean(ss[,"a"]) > -0.8, msg="did not shift means(a) towards narrow part." )
		})	
}
#twUtestF(getSubSpaces)

