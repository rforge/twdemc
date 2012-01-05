#generating a twDEMC result with den2dCor

(.expTheta <- c(a=0,b=0) )
(.expCovTheta <- diag(c(a=2,b=2)) )		
.nPop=2
Zinit <- initZtwDEMCNormal( .expTheta, .expCovTheta, nChainPop=4, nPop=.nPop)
#mtrace(twDEMCBlockInt)

argsFLogDen = list()
do.call( den2dCor, c(list(theta=Zinit[,1,1]),argsFLogDen))

den2dCorTwDEMCPops <- thin( twDEMCBlock(Zinit, nGen=1000, dInfos=list(list(fLogDen=den2dCor)), nPop=.nPop, debugSequential=TRUE ), start=300)
den2dCorTwDEMC <- concatPops(den2dCorTwDEMCPops)
#den2dCorTwDEMC <- twDEMCBatch(den2dCorTwDEMC, nGen=1500)
#den2dCorTwDEMC3 <- twDEMCBatch(den2dCorTwDEMC, nGen=1000+6*500)	# compare to divideTwDEMC
#str(den2dCorTwDEMC)
getNGen(den2dCorTwDEMCPops)
getSpacesPop(den2dCorTwDEMCPops)

.tmp.f <- function(){
	#den2dCorTwDEMC <- concatPops(den2dCorTwDEMC)
	rescoda <- as.mcmc.list(den2dCorTwDEMC)
	#rescoda <- as.mcmc.list(den2dCorTwDEMC3)
	plot(rescoda, smooth=FALSE)
	ss <- stackChains(den2dCorTwDEMC)
	#ss <- stackChains(den2dCorTwDEMC3)
	plot(density(ss[,"a"]))
	plot( b ~ a, as.data.frame(ss) ); 
	plot( b ~ a, as.data.frame(ss), xlim=c(-0.5,2), ylim=c(-20,40) ); 
	points(0.8, 0, col="red" )	# theoretical maximum
}

#-----  infer subspaces on a subsample
tmp <- squeeze(den2dCorTwDEMCPops, length.out= 256 %/% getNChainsPop(den2dCorTwDEMCPops) ) # thin to 256 samples per space
ss0 <- stackChainsPop(concatPops(tmp))
#plot( b ~ a, as.data.frame(ss0[,,1]), xlim=c(-0.5,2), ylim=c(-20,40) )

nBlock <- attr(ss0,"nBlock")
minPSub <- 0.1
den2dCorSubSpaces <- lapply( 1:dim(ss0)[3], function(iPop){
		samplePop <- ss0[,,iPop]
		#mtrace(getSubSpaces)
		#mtrace(findSplit)
		getSubSpaces(samplePop[,-(1:nBlock)], minPSub=minPSub, isBreakEarly=FALSE, argsFSplit=list())	# here omit the logDensity column
	})

#------- single step MC runs using subspaces
den2dCorTwDEMCSpaces <-  divideTwDEMCPops(den2dCorTwDEMCPops, den2dCorSubSpaces )
getSpacesPop(den2dCorTwDEMCSpaces)

# ------ iterative MC runs using subspaces
# XXTODO

# save only the 1000 generations run, it can be easily extended
den2dCorEx <- list(
	den2dCor = den2dCor				##<< the density function
	,mcBulk = den2dCorTwDEMCPops	##<< DEMC run (class twDEMCPops) without subspaces
	,subspaces0 = den2dCorSubSpaces	##<< subspaces inferred on a subsample of mcBulk
	,mcSubspaces0 = den2dCorTwDEMCSpaces ##<< DEMC run using a single step with subspaces 
	)
save(den2dCorEx  , file="data/den2dCorEx.RData")


.tmp.f <- function(){
	aSample <- stackChainsPop(den2dCorTwDEMC)
	den2dCorDivideTwDEMC <- divideTwDEMCBatch( aSample, nGen=8*512, fLogDen=den2dCor)
	.tmp.f <- function(){
		rescoda <- as.mcmc.list( lapply(twListArrDim(den2dCorDivideTwDEMC$sample),mcmc)  )
		plot(rescoda, smooth=FALSE)
		ss <- stackChains(den2dCorTwDEMC)
		ss2 <- t(stackChains(den2dCorDivideTwDEMC$sample))
		plot(density(ss[,"a"]))
		lines(density(ss2[,"a"]), col="green")
		plot( b ~ a, as.data.frame(ss2), xlim=c(-0.5,2), ylim=c(-20,40) ); 
		points(0.8, 0, col="red" )	# theoretical maximum
	}
	save(den2dCorDivideTwDEMC ,file="data/den2dCorDivideTwDEMC.RData")
}

