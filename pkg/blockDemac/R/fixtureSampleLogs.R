.fixtureSampleLogs <- function( 
    sampleDimensions=.fixtureSampleDimensions()
    ,x0 = c(a=1, b=7, c=9)                # state at first step
    ,logDensityComponents0 = c(obs1=-10, obs2=-5, parms=-2) # logDensity at first step
    ,countAcceptedInInterval0 = structure(rep(1L,getNBlock(sampleDimensions)),names=getBlockNames(sampleDimensions))
    ,proportionAcceptedInInterval0 = structure(rep(0.25,getNBlock(sampleDimensions)),names=getBlockNames(sampleDimensions))
    ,thin=4L
    ,nInitialSample=0L
){
    nPopulation <- getNPopulation(sampleDimensions)
    nSample <- getNSamplePopulations(sampleDimensions)
    nChainInPopulation <- getNChainInPopulation(sampleDimensions)
    # if x is a vector, generate x0 for different chains by multiplying by chainInPopulation index
    if(!is.matrix(x0)) x0 <- abind(lapply( 1:nPopulation, function(iPop){
                    abind( lapply( 1:nChainInPopulation, function(iChain){
                            x0 * iChain                                
                        }), along=2)
                }))
    if(!is.matrix(logDensityComponents0))
            logDensityComponents0 <- abind(lapply( 1:nPopulation, function(iPop){
                        abind(lapply( 1:nChainInPopulation, function(iChain){
                                logDensityComponents0 * iChain                                
                            }),along=2)
                    }),along=2)
    # for proportionAcceptedInIntervalChains, multiply by chainInPop index modulo thin
    if(!is.matrix(countAcceptedInInterval0)) countAcceptedInInterval0 <- abind(lapply( 1:nPopulation, function(iPop){
                    abind( lapply( 1:nChainInPopulation, function(iChain){
                                #TODO replace countAccepted by pAccept to remove dependency on thin 
                            countAcceptedInInterval0 *  (((iChain-1) %% thin)+1)                                
                        }), along=2)
                }))
    if(!is.matrix(proportionAcceptedInInterval0)) proportionAcceptedInInterval0 <- abind(lapply( 1:nPopulation, function(iPop){
                    abind( lapply( 1:nChainInPopulation, function(iChain){
                                #TODO replace countAccepted by pAccept to remove dependency on thin 
                                proportionAcceptedInInterval0 *  (((iChain-1) %% thin)+1)                                
                            }), along=2)
                }))
    if( nrow(x0) != getNParameter(sampleDimensions) ) stop("mismatch between x0 and number of parameters in sampleDimensions.")        
    if( !all(rownames(x0) == getParameterNames(sampleDimensions)) ) stop("mismatch between names in x0 and sampleDimensions.")     
    if( nrow(logDensityComponents0) != getNLogDensityComponent(sampleDimensions)) stop("mismatch between logDensityComponents0 and number of logDensity components in sampleDimensions.")
    if( !all(rownames(logDensityComponents0) == getLogDensityComponentNames(sampleDimensions)) ) stop("mismatch between names in logDensityComponents0 and sampleDimensions.")     
    nChain <- ncol(x0)
    assert_that( nChain == nChainInPopulation * nPopulation )
    assert_that( nChain == ncol(logDensityComponents0) )
    assert_that( nChain == ncol(countAcceptedInInterval0) )
    assert_that( nChain == ncol(proportionAcceptedInInterval0) )
    #
    populationLogs <- lapply(1:nPopulation, function(iPop){
            # multiply by step
            parametersChains <- abind( lapply(1:nChainInPopulation, function(iChain){
                        x0Chain <- x0[,(iPop-1)*nChainInPopulation + iChain]
                        abind( lapply( 1:nSample[iPop], function(i){ i*x0Chain }), along=2)
                    }), rev.along=0)
            logDensityComponentsChains <- 
                    abind( lapply(1:nChainInPopulation, function(iChain){
                                logDensityComponents0Chain <- logDensityComponents0[,(iPop-1)*nChainInPopulation + iChain]
                                abind( lapply( 1:nSample[iPop], function(i){ i*logDensityComponents0Chain }),along=2)
                            }), rev.along=0)
            # for countAcceptedInIntervalChains do not multiply by step, but all steps the same
            countAcceptedInIntervalChains <- abind( lapply(1:nChainInPopulation, function(iChain){
                        countAcceptedInInterval0Chain <- countAcceptedInInterval0[,(iPop-1)*nChainInPopulation + iChain]
                        abind( lapply( 1:nSample[iPop], function(i){ 1*countAcceptedInInterval0Chain }), along=2)
                    }), rev.along=0)
            # for proportionAcceptedInIntervalChains do not multiply by step, but all steps the same
            proportionAcceptedInIntervalChains <- abind( lapply(1:nChainInPopulation, function(iChain){
                        proportionAcceptedInInterval0Chain <- proportionAcceptedInInterval0[,(iPop-1)*nChainInPopulation + iChain]
                        abind( lapply( 1:nSample[iPop], function(i){ 1*proportionAcceptedInInterval0Chain }), along=2)
                    }), rev.along=0)
            #
            stepTemperatures <- matrix(nrow=nrow(logDensityComponents0), ncol=nSample[iPop], dimnames=list(rownames(logDensityComponents0),NULL))
            stepTemperatures[] <- 1
            # if nSample==0, then row index is 1:0, need to subset to include zero rows
            isAtLeastOneSample <- (nSample[iPop] != 0)
            if( ncol(parametersChains) < nInitialSample ) stopDemac("requested too many fixture initial parameters.")
            # reuse parameters for initial sample (used for initializing the jumps)
            parametersChainsWithInitial <- if( nInitialSample == 0L) parametersChains else
                abind(parametersChains[,1:nInitialSample, ,drop=FALSE], parametersChains, along=2L)            
            sampleLog <- new("SampleLogPopulation"
                , parameters=parametersChainsWithInitial[,isAtLeastOneSample, ,drop=FALSE]
                , logDensityComponents=logDensityComponentsChains[,isAtLeastOneSample, ,drop=FALSE]
                , stepTemperatures=stepTemperatures
                , nInitialSample=nInitialSample
            )
        })
    sampleLogs <- new("SampleLogs", sampleDimensions=sampleDimensions, populationLogs=populationLogs )
    list(
        x0 = x0[,1]
        ,x0Chains = x0
        ,logDensityComponents0 = logDensityComponents0[,1]
        ,logDensityComponents0Chains = logDensityComponents0
        # TODO remove countAcceptedInInterval
        ,countAcceptedInInterval0 = countAcceptedInInterval0[,1]
        ,countAcceptedInInterval0Chains = countAcceptedInInterval0
        ,proportionAcceptedInInterval0 = proportionAcceptedInInterval0[,1]
        ,proportionAcceptedInInterval0Chains = proportionAcceptedInInterval0
        ,sampleDimensions=sampleDimensions
        ,sampleLogs = sampleLogs
    )
}
attr(.fixtureSampleLogs,"ex") <- function(){
    sampleLogs <- .fixtureSampleLogs()
    getParametersForPopulation(sampleLogs$sampleLog,1L)     # note 26 steps, x0 * iChainInPop * iStep
    getLogDensityComponentsForPopulation(sampleLogs$sampleLog,2L)  # note 30 steps, logDensity0 * iChainInPop * iStep
}

.fixtureSampleLogsMetr <- function(
    fx=.fixtureChainSampler()
    ,sampleDimensions = .fixtureSampleDimensions(blockSpecifications=fx$blocks)
    ,x0 = fx$x0
    ,logDensityComponents0 = fx$logDensityComponents0
    # TODO remove countAcceptedInInterval
    ,countAcceptedInInterval0 = c(add7=1,met1=2,met2=3)        # number of accepted steps at first sample, ie. thinning interval
    ,...
){
    sampleLogs <- .fixtureSampleLogs(
        sampleDimensions = sampleDimensions
        ,x0 = x0
        ,logDensityComponents0 = logDensityComponents0
        ,countAcceptedInInterval0 = countAcceptedInInterval0
        ,...
    )
}
attr(.fixtureSampleLogsMetr,"ex") <- function(){
    sampleLogs <- .fixtureSampleLogsMetr()
    getParametersForPopulation(sampleLogs$sampleLog,1L)     # note 26 steps, x0 * iChainInPop * iStep
    getLogDensityComponentsForPopulation(sampleLogs$sampleLog,2L)  # note 30 steps, logDensity0 * iChainInPop * iStep
}
