.setUp <- function(){
	#attach(twLinreg1)
	#attach(twdemcEx1)
}

.tearDown <- function(){
	#detach( twdemcEx1 )
	#detach( twLinreg1 )
}


test.divideTwDEMCStep <- function(){
    data(den2dCorEx)
    .nGen <- 256
    aTwDEMC <- den2dCorEx$mcSubspaces0
    m0 <- calcM0twDEMC( getNParms(aTwDEMC), getNChainsPop(aTwDEMC) )
    .nRowsMin <- 2*(m0+1)
	res <- divideTwDEMCStep(aTwDEMC, nGen=.nGen, dInfos=list(list(fLogDen=den2dCor))
            , 	nRowsMin = .nRowsMin		# minimum number of rows in population, so that enough samples to split into 2 subs
            ,   nRowsMax = .nRowsMin		# minimum number of rows in population, so that enough samples to split into 2 subs
            ,   minPSub = 0.1
            ,  TEnd = getCurrentTemp(aTwDEMC)
        , debugSequential=TRUE )
	mc1 <- stackPopsInSpace(res$resTwDEMC) 
	.tmp.f <- function(){
        getNSamples(res$resTwDEMC)
        getNSpaces(res$resTwDEMC)
        
        #windows(record=TRUE)
        plot( as.mcmc.list(stackPopsInSpace(aTwDEMC)), smooth=FALSE)
        #mc1 <- stackPopsInSpace(res$resTwDEMC, mergeMethod="stack")
        mc1 <- stackPopsInSpace(res$resTwDEMC, mergeMethod="random")
        plot( as.mcmc.list(mc1), smooth=FALSE) # note how the distribution shifted
        #plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=4)), smooth=FALSE)
        #plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=12)), smooth=FALSE)
        
        ss <- stackChains(concatPops(mc1))			# the new sample
        ss0 <- stackChains(concatPops(squeeze(stackPopsInSpace(aTwDEMC),length.out=getNSamples(mc1)[1]))) # initial sample of same length
        plot( b ~ a, as.data.frame(ss0), xlim=c(-0.5,2), ylim=c(-20,40) )
        plot( b ~ a, as.data.frame(ss), xlim=c(-0.5,2), ylim=c(-20,40) ) # note that more samples are in the region of interest
        points(0.8, 0, col="red" )	# theoretical maximum
        
        
	}
	checkEquals( getNSpaces(den2dCorEx$mcSubspaces0), getNPops(mc1), msg="wrong number of spaces" )
	
	iPop=1
	sapply(1:getNPops(mc1), function(iPop){
		pop <- mc1$pops[[iPop]]
		checkTrue( nrow(pop$parms) >  256/4-10, msg=" wrong number of generations in subPopulation") # number of generations 
		ss <- stackChains(subPops(mc1, iPops=iPop ))
		imax <- which.max(ss[,1])
		thetaHat <- ss[imax,]			
		checkInterval(thetaHat["a"], -0.8, +1.5,msg=paste("wrong thetaHat[a]=",thetaHat["a"],sep="") )
		checkInterval(thetaHat["b"], -20, +20,msg=" wrong thetaHat[b]")
		checkTrue( mean(ss[,"a"]) > -0.8, msg="did not shift means(a) towards narrow part." )
	})	
}
#twUtestF(getSubSpaces)

# see example of divideTwDEMCSACon


