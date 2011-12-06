#generating a twDEMC result with den2dCor

(.expTheta <- c(a=0,b=0) )
(.expCovTheta <- diag(c(a=2,b=2)) )		
.nPop=2
Zinit <- initZtwDEMCNormal( .expTheta, .expCovTheta, nChainPop=4, nPop=.nPop)
#mtrace(twDEMCInt)

argsFLogDen = list()
do.call( den2dCor, c(list(theta=Zinit[,1,1]),argsFLogDen))

den2dCorTwDEMC <- concatPops(twDEMCBlock(Zinit, nGen=1000, dInfos=list(list(fLogDen=den2dCor)), nPop=.nPop, debugSequential=TRUE ))
den2dCorTwDEMC <- twDEMCBatch(den2dCorTwDEMC, nGen=1500)
den2dCorTwDEMC3 <- twDEMCBatch(den2dCorTwDEMC, nGen=1000+6*500)	# compare to divideTwDEMC
str(den2dCorTwDEMC)

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

# save only the 1000 generations run, it can be easily extended
save(den2dCorTwDEMC ,file="data/den2dCorTwDEMC.RData")


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

