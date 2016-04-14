.prof.proposeSampleAndLogRange <- function(){
    #library(profr)
    #library(lineprof)    
    object0 <- object
    #p1 <- profr({ for( i in 1:2 ){
    p1 <- profr({ for( i in 1:100 ){
        object <- object0
        iChainsStep <- getIChainsForPopulations(getSampleDimensions(object@sampleLogs), iPopulationsStep)
        pAcceptPops  <- computePopulationAcceptanceRates(object@acceptanceTracker)
        iPopulationChains <- getIPopulationsForChains(getSampleDimensions(object@sampleLogs), iChainsStep)
        pAcceptChains <- pAcceptPops[iPopulationChains]
        # jumps relate already only to the populations taking part in this range (iPopulationsStep)
        jumps <- proposeJumps( object@jumpProposer
                , nGeneration=nThinningIntervalStep*getThin(object@chainSampler)
                , sampleLogs=object@sampleLogs
                , iPopulationsStep = iPopulationsStep
                , iCurrentSamples = object@nSampleBeforeSampling+iSampleInBatch-1  # current state before the current index
                , acceptanceRates = pAcceptPops
        )
        # samplesChains also takes information of the participating populations
#        samplesChains <- sampleRangeAllChains(object
#                , chainStates = object@chainStates[iChainsStep]
#                , intervalInfoChains =list(
#                        iChainsStep=iChainsStep
#                        ,iSample0=iSampleInBatch-1
#                        ,nSample=nThinningIntervalStep
#                        ,step=jumps$jump
#                        , rExtra=jumps$logSnookerDensityMultiplier
#                        , pAccept=pAcceptChains
#                ), mcCommon=.getMcCommon(object))
        # samplesChains holds only the participating chains
        object@chainStates[iChainsStep] <- lapply( samplesChains, "[[", "chainState" )
        sampleLogChains <- lapply( samplesChains, "[[", "sampleLog" )
        object@sampleLogs <- recordSampleLogChains(object@sampleLogs, sampleLogChains
                , iPopulations=iPopulationsStep, iSample0=object@nSampleBeforeSampling+iSampleInBatch-1 )
        proportionAcceptedInInterval <- abind( lapply(sampleLogChains, getProportionAcceptedInInterval), rev.along=0)
        object@acceptanceTracker <- recordAcceptanceRate(object@acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval, iChains=iChainsStep)
        
    #} })
    #} }, interval = 0.0005)
    } }, interval = 0.01)
plot(p1)
}
# mostly .sampleJumpsFromHistory



# mostly .sampleJumpsFromHistory



.prof.pcolsEqual <- function(){
    # in .sampleStates
    library(rbenchmark)
    benchmark(
        iSame <- twWhichColsEqual( adrop(zx[,,2 ,drop=FALSE],3), adrop(zx[,,1 ,drop=FALSE],3) )
        ,iSame <- twWhichColsEqual( zx[,,2], zx[,,1] ) 
        ,iSame <- whichColsEqualSumHeuristics( zx[,,2], zx[,,1] ) 
        ,replications=1000 	
    )
    # speedup of factor 2 with each optimization
}

.prof.abind1 <- function(){
    # in .sampleStates
    benchmark(
            zx0 <- abind( zParms, Xs, along=3)
            ,zx <- array( c(zParms, Xs), dim=c(dim(zParms)[1:2], 4), dimnames=c(dimnames(zParms)[1], NULL,NULL) )
        ,replications=1000 	
    )
    #speedup of 8
}

.prof.abind2 <- function(){
    # in .sampleJumpsFromHistory
    benchmark(
            zx0 <- abind( zxPops, along= 2 )	# combine chains (and steps within each chain) of all populations
            ,zx <- abind3DAlong2(zxPops)
            ,replications=1000 	
    )
    #speedup of 2.8
}




