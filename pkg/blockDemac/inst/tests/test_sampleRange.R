#require(testthat)
context("sampleRange")


fUpdateBlockAddTerm <- function(
    x, iParametersToUpdate
    , lowerParBounds, upperParBounds, intermediates
    , term
){
    list(
        isUpdated = TRUE 
        ,xC = x[ iParametersToUpdate ] + term
        ,intermediate = NULL    # actually dont return NULL, better list(), here testing handled
    )
}

fLogDenConst <- function(
    x                           ##<< parameter vector
    ,intermediates              ##<< intermediate results
    ,logDensityComponents       ##<< numeric vector (nLDensityComponents): 
        ## already known logDensitys. Only calculate the positions that are NA
    ,...
    ,acceptedLogDensityComponents
    ,logEnv=NULL        # usually in caller: new.env(parent = emptyenv())
){
    # if first parameter is 1 then accept, else not
    if( is.environment(logEnv) ){
        logEnv$log <- paste0(logEnv$log,",x=",catNamedVector(x,3)," logDenC=",catNamedVector(logDensityComponents))
    }
    ret <- if (x[1]==1) {
            acceptedLogDensityComponents    # accepted  
        } else {
            Inf*acceptedLogDensityComponents    # not accepted
        }
    ret
}

fLogDen1D <- function(
        x                           ##<< parameter vector
        ,intermediates              ##<< intermediate results
        ,logDensityComponents       ##<< numeric vector (nLDensityComponents): 
            ## already known logDensitys. Only calculate the positions that are NA
    ,...
    ,obs=5
    ,sdObs=1
){
    # sample x[1] around obs=5 
    c(obs1=-1/2* ((x[1]-obs)/sdObs)^2)
}

x0 <- c(a=0, isM=1, m=0, isC=1, c=0)
chainState <- new("ChainState", 
    parameters=x0, 
    logDensityComponents=c(m1=-100, m2=-100, parms=-200)
)

thin <- 4L
nInterval <- 3L
nGen <- thin*nInterval

stepSingle <- c(a=NA, isM=0, m=7, isC=0, c=9)   # isX stays at initial 1
stepInfoRange <- new("StepInfoRangeImpl", step=matrix(
        rep(stepSingle, nGen), nrow=length(x0), ncol=nGen, dimnames=list(names(x0),NULL))
    ,nLogDensityComponents=length(getLogDensityComponents(chainState)))

blocks <- list(
    add7 = blockSpec("a",, new("FunctionBasedBlockUpdater", 
            fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=7))),
    # important that isM is the first parameter used
    met1 = blockSpec(c("isM","m"),c("isM","m","c"),new("MetropolisBlockUpdater",
            fLogDensity=fLogDenConst,
            argsFLogDensity=list(acceptedLogDensityComponents = c(m1=-10)),
            logDensityComponentNames = c("m1"))), 
    # important that isC is the first parameter used
    met2 = blockSpec(c("isC","c"),c("isC","c","m"),new("MetropolisBlockUpdater",
            fLogDensity=fLogDenConst,
            argsFLogDensity=list(acceptedLogDensityComponents = c(m2=-10,parms=-2)),
            logDensityComponentNames = c("m2","parms"))) 
)

#blockDimensions <- newBlockDimensions(blocks, names(x0))
blockUpdaters <- newBlockUpdaters(blocks, list(), names(x0))

chainSampler <- new("ChainSamplerImpl", blockUpdaters=blockUpdaters, thin=thin)
nSamplePopulations <- c(16L,32L)
nChainInPopulation <- 2L
sampleDimensions <- initializeSampleDimensionsImpl( new("SampleDimensionsImpl")
    , blockDimensions=getBlockDimensions(blockUpdaters)
    , nSamplePopulations=nSamplePopulations, nChainInPopulation=nChainInPopulation)


mcCommon <- list(
    sampleDimensions = sampleDimensions
    ,chainSampler = chainSampler
    ,subSpacePopulations = rep(list( initializeIndexedBoundsForParameterNames(new("SubSpace"),names(x0))), 4)
    ,temperatures = rep(list(pop1=stepInfoRange@temperature),4)
)

populationSampler <- new("PopulationSampler")        

test_that(".sampleRangeOneChain",{
        res <- .sampleRangeOneChain(chainState, iPopulation=1
            , iSample0=0L, nSample=nInterval
            , stepInfoRange, mcCommon)
        updatedChainState <- res$chainState
        sampleLog <- res$sampleLog
        parameters <- getParameters(sampleLog)
        logDensityComponents <- getLogDensityComponents(sampleLog)
        proportionAcceptedInInterval <- getProportionAcceptedInInterval(sampleLog)
        #
        expect_that( parameters[,3], equals(x0 +3*c(a=7*thin, isM=0, m=7*thin, isC=0, c=9*thin)))
        expect_that( proportionAcceptedInInterval[,3], is_equivalent_to(rep(thin/thin,3)))
        expect_that( logDensityComponents[,3], equals(c(m1=-10,m2=-10,parms=-2)))
    })

test_that(".sampleRangeAllChains",{
        .tmp.f.depr <- function(){
            sfInit(parallel=TRUE, cpus=2)
            sfLibrary(blockDemac)
            sfExport("mcCommon", ".sampleRangeOneChain", local=FALSE)
        }
        .tmp.f <- function(){
            cl <- makeCluster(2L)
            clusterExport(cl, c("mcCommon", ".sampleRangeOneChain"), envir=environment())
            populationSampler@cluster <- cl
            # when omitting loading library, gives an error (check if really done on cluster)
            clusterEvalQ(cl, library(blockDemac) ) 
            #
            .tmp.f <- function(){
                clusterEvalQ(cl, helloFromBlockDemac() )
                dumpedStop <- dumpOnError(stop)
                lapply( 1, dumpedStop("dumpedStop"))
                load("last.dump.rda")
                debugger()
            }
            #stopCluster(cl)
        }
        
        nChain <- getNChain(sampleDimensions)
        chainStates <- rep(list(chainState),nChain)
        iPopulationChains <- c(1,1,2,2)
        # step for each chain is multiplied by chainNumber (1:nChain)
        step <- abind( lapply(1:nChain, function(i){ stepInfoRange@step*i }), rev.along=0)
        rExtra <- abind( rep(list(stepInfoRange@rExtra),nChain), rev.along=0)
        pAccept <- rep(0.25, nChain)
        res <- sampleRangeAllChains(populationSampler, chainStates, intervalInfoChains =list(
                #iPopulationChains=iPopulationChains
                iChainsStep=1:4
                , iSample0=0L, nSample=nInterval
                ,step=step, rExtra=rExtra, pAccept=pAccept
            ), mcCommon=mcCommon)
        expect_that( length(res), equals(nChain) )
        #
        updatedChainState <- res[[1]]$chainState
        sampleLog <- res[[1]]$sampleLog
        parameters <- getParameters(sampleLog)
        logDensityComponents <- getLogDensityComponents(sampleLog)
        proportionAcceptedInInterval <- getProportionAcceptedInInterval(sampleLog)
        expect_that( parameters[,3], equals(x0 +3*c(a=7*thin, isM=0, m=7*thin, isC=0, c=9*thin)))
        expect_that( proportionAcceptedInInterval[,3], is_equivalent_to(rep(thin/thin,3)))
        expect_that( logDensityComponents[,3], equals(c(m1=-10,m2=-10,parms=-2)))
        #
        updatedChainState <- res[[3]]$chainState
        sampleLog3 <- res[[3]]$sampleLog
        parameters <- getParameters(sampleLog3)
        logDensityComponents <- getLogDensityComponents(sampleLog3)
        proportionAcceptedInInterval <- getProportionAcceptedInInterval(sampleLog3)
        expect_that( parameters[-1,3], equals(x0[-1] +3*3*c(isM=0, m=7*thin, isC=0, c=9*thin)))
        expect_that( proportionAcceptedInInterval[,3], is_equivalent_to(rep(thin/thin,3)))
        expect_that( logDensityComponents[,3], equals(c(m1=-10,m2=-10,parms=-2)))
    })

test_that(".sampleRangeAllChains zero logDen",{
        x0 <- c(a=0)
        chainState <- new("ChainState", 
            parameters=x0, 
            logDensityComponents=numeric(0)
        )
        stepSingle <- c(a=0)  
        stepInfoRange <- new("StepInfoRangeImpl", step=matrix(
                rep(stepSingle, nGen), nrow=length(x0), ncol=nGen, dimnames=list(names(x0),NULL))
            ,nLogDensityComponents=length(getLogDensityComponents(chainState)))
        blocks <- list(
            add7 = blockSpec("a",, new("FunctionBasedBlockUpdater", 
                    fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=7)))
        )
        blockUpdaters <- newBlockUpdaters(blocks, list(), names(x0))
        chainSampler <- new("ChainSamplerImpl", blockUpdaters=blockUpdaters, thin=thin)
        nSamplePopulations <- c(2L,3L)
        sampleDimensions <- initializeSampleDimensionsImpl( new("SampleDimensionsImpl")
            , blockDimensions=getBlockDimensions(blockUpdaters)
            , nSamplePopulations=nSamplePopulations, nChainInPopulation=nChainInPopulation)
        mcCommon <- list(
            sampleDimensions = sampleDimensions
            ,chainSampler = chainSampler
            ,subSpacePopulations = rep(list( initializeIndexedBoundsForParameterNames(new("SubSpace"),names(x0))), 4)
            ,temperatures = rep(list(pop1=stepInfoRange@temperature),4)
        )
        #
        nChain <- getNChain(sampleDimensions)
        chainStates <- rep(list(chainState),nChain)
        iPopulationChains <- c(1,1,2,2)
        # step for each chain is multiplied by chainNumber (1:nChain)
        step <- abind( lapply(1:nChain, function(i){ stepInfoRange@step*i }), rev.along=0)
        rExtra <- abind( rep(list(stepInfoRange@rExtra),nChain), rev.along=0)
        pAccept <- rep(0.25, nChain)
        res <- sampleRangeAllChains(populationSampler, chainStates, intervalInfoChains =list(
                #iPopulationChains=iPopulationChains
                iChainsStep=1:4
                , iSample0=0L, nSample=nInterval
                ,step=step, rExtra=rExtra, pAccept=pAccept
            ), mcCommon=mcCommon)
        expect_that( length(res), equals(nChain) )
        #
        updatedChainState <- res[[3]]$chainState
        sampleLog <- res[[3]]$sampleLog
        parameters <- getParameters(sampleLog)
        logDensityComponents <- getLogDensityComponents(sampleLog)
        proportionAcceptedInInterval <- getProportionAcceptedInInterval(sampleLog)
        expect_that( parameters[,3], equals(x0 +3*c(a=7*thin)))
        expect_that( proportionAcceptedInInterval[,3], is_equivalent_to(rep(thin/thin,1)))
        expect_that( logDensityComponents[,3], equals(numeric(0)))
    })

test_that(".sampleRangeAllChains 1D logDen",{
        x0 <- c(a=-2)
        chainState <- new("ChainState", 
            parameters=x0, 
            logDensityComponents=c(obs1=NA_real_)
        )
        thin <- 2L
        nInterval <- 4L #24L
        nGen <- thin*nInterval
        stepInfoRange <- new("StepInfoRangeImpl", step=matrix(
                rnorm(nGen, 0, 1.2)
                , nrow=length(x0), ncol=nGen, dimnames=list(names(x0),NULL))
            ,nLogDensityComponents=length(getLogDensityComponents(chainState)))
        blocks <- list(
            met1 = blockSpec(,,new("MetropolisBlockUpdater",
                    fLogDensity=fLogDen1D,
                    argsFLogDensity=list(),
                    logDensityComponentNames = c("obs1"))) 
        )
        blockUpdaters <- newBlockUpdaters(blocks, list(), names(x0))
        chainSampler <- new("ChainSamplerImpl", blockUpdaters=blockUpdaters, thin=thin)
        nSamplePopulations <- c(2L,nInterval)
        sampleDimensions <- initializeSampleDimensionsImpl( new("SampleDimensionsImpl")
            , blockDimensions=getBlockDimensions(blockUpdaters)
            , nSamplePopulations=nSamplePopulations, nChainInPopulation=nChainInPopulation)
        mcCommon <- list(
            sampleDimensions = sampleDimensions
            ,chainSampler = chainSampler
            ,subSpacePopulations = rep(list( initializeIndexedBoundsForParameterNames(new("SubSpace"),names(x0))), 4)
            ,temperatures = rep(list(pop1=stepInfoRange@temperature),4)
        )
        #
        nChain <- getNChain(sampleDimensions)
        chainStates <- rep(list(chainState),nChain)
        iPopulationChains <- c(1,1,2,2)
        # step for each chain is multiplied by chainNumber (1:nChain)
        step <- abind( lapply(1:nChain, function(i){ stepInfoRange@step*i }), rev.along=0)
        rExtra <- abind( rep(list(stepInfoRange@rExtra),nChain), rev.along=0)
        pAccept <- rep(0.25, nChain)
        res <- sampleRangeAllChains(populationSampler, chainStates, intervalInfoChains =list(
                #iPopulationChains=iPopulationChains
                iChainsStep=1:4
                , iSample0=0L, nSample=nInterval
                ,step=step, rExtra=rExtra, pAccept=pAccept
            ), mcCommon=mcCommon)
        expect_that( length(res), equals(nChain) )
        #
        iChain <- 3L
        iPop <- getIPopulationsForChains(sampleDimensions,iChain)
        updatedChainState <- res[[3]]$chainState
        sampleLog <- res[[3]]$sampleLog
        parameters <- getParameters(sampleLog)
        logDensityComponents <- getLogDensityComponents(sampleLog)
        proportionAcceptedInInterval <- getProportionAcceptedInInterval(sampleLog)
        expect_that( dim(parameters), equals(c(1, nSamplePopulations[iPop])))
        expect_that( dim(logDensityComponents), equals(c(1,nSamplePopulations[iPop])))
        .tmp.f <- function(){
            plot(1:nInterval, as.vector(parameters), type="l" )
        }
    })


# commented because call to sfLibrary(blockDemac)
#test_that("StepInfoRange in remote process (learning test)",{
#        x0 <- c(a=0, isM=1, m=0, isC=1, c=0)
#        chainState <- new("ChainState", 
#            parameters=x0, 
#            logDensityComponents=c(m1=-100, m2=-100, parms=-200)
#        )
#        .nGen <- 2
#        stepSingle <- c(a=NA, isM=0, m=7, isC=0, c=9)   # isX stays at initial 1
#        stepInfoRange <- new("StepInfoRangeImpl", step=matrix(
#                rep(stepSingle, .nGen), byrow=TRUE, nrow=.nGen, ncol=length(x0), dimnames=list(NULL,names(x0)))
#            ,nLogDensityComponents=length(getLogDensityComponents(chainState)))
#        #        
#        sfInit(parallel=TRUE, cpus=2)
#        sfLibrary(blockDemac)
#        fRemote <- function(stepInfoRange){ getStep(stepInfoRange) }
#        res <- sfClusterCall(fRemote, stepInfoRange=stepInfoRange)
#        tmp <- lapply(res, function(resNode){
#            expect_that(resNode, equals(stepSingle))                
#            })
#    })



