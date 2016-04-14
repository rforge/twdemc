#require(testthat)
context(".initChainStates")

thin <- 4L
fx <- .fixtureSampleLogs(
    sampleDimensions = .fixtureSampleDimensions(
        nSample = c(24L,32L)
        ,nChainInPopulation = 4L
    )
)
sampleLogs <- fx$sampleLogs
sDim <- getSampleDimensions(sampleLogs)
# sampleLog will be x0*iChain*iStep
# countAcceptedInInterval=1L for all 


# block over entire parameter vector
blockIndices <- new("BlockIndices"
    , iParametersToUpdate=1:getNParameter(sDim)
    , iParametersUsed=1:getNParameter(sDim)
    , intermediateIdsUsed="nonExisting" 
)

chainStates0 <-  lapply(1:getNChain(sDim), function(iChain){ 
        new("ChainState", parameters=fx$x0Chains[,iChain]/2, logDensityComponents=fx$logDensityComponents0Chains[,iChain]/2)
    })



test_that(".initChainStatesOneChain for first row",{
        logDensityComponents0 <- fx$logDensityComponents0                
        chainStates <- .initChainStates(sampleLogs,iSamplePopulations=rep(1L,getNPopulation(sDim)), chainStates=chainStates0)
        expect_that( length(chainStates), equals(getNChain(sDim)))
        iChain <- 1L
        chainState <- chainStates[[iChain]]
        expect_that( getChainStateParameters(chainState), equals(fx$x0Chains[,1]*iChain))
        expect_that( getLogDensityComponents(chainState), equals(logDensityComponents0*iChain))
        expect_that( isChangedByLastUpdate(chainState,blockIndices), equals(rep(FALSE,nrow(fx$x0Chains))) )
        #expect_true( is.null(getBlockIntermediate(chainState, blockIndices)) )
        iChain <- 2L    # 3L is modulo nPop equal to 1
        chainState <- chainStates[[iChain]]
        expect_that( getChainStateParameters(chainState), equals(fx$x0Chains[,1]*iChain))
        expect_that( getLogDensityComponents(chainState), equals(logDensityComponents0*iChain))
        expect_that( isChangedByLastUpdate(chainState,blockIndices), equals(rep(FALSE,nrow(fx$x0Chains))) )
        #expect_true( is.null(getBlockIntermediate(chainState, blockIndices)) )
    })

test_that(".initChainStatesOneChain for third row",{
        logDensityComponents0 <- fx$logDensityComponents0                
        iSample <- 3L
        chainStates <- .initChainStates(sampleLogs, iSamplePopulations=rep(iSample, getNPopulation(sDim)), chainStates=chainStates0)
        expect_that( length(chainStates), equals(getNChain(sDim)))
        iChain <- 1L
        chainState <- chainStates[[iChain]]
        expect_that( getChainStateParameters(chainState), equals(fx$x0Chains[,1]*iChain*iSample))
        expect_that( getLogDensityComponents(chainState), equals(logDensityComponents0*iChain*iSample))
        expect_that( isChangedByLastUpdate(chainState,blockIndices), equals(rep(FALSE,nrow(fx$x0Chains))) )
        #expect_true( is.null(getBlockIntermediate(chainState, blockIndices)) )
        iChain <- 2L    # 3L is modulo nPop equal to 1
        chainState <- chainStates[[iChain]]
        expect_that( getChainStateParameters(chainState), equals(fx$x0Chains[,1]*iChain*iSample))
        expect_that( getLogDensityComponents(chainState), equals(logDensityComponents0*iChain*iSample))
        expect_that( isChangedByLastUpdate(chainState,blockIndices), equals(rep(FALSE,nrow(fx$x0Chains))) )
        #expect_true( is.null(getBlockIntermediate(chainState, blockIndices)) )
    })

test_that(".initChainStatesOneChain default to last row of the sample",{
        logDensityComponents0 <- fx$logDensityComponents0                
        iSample <- min(getNSamplePopulations(sampleLogs))
        chainStates <- .initChainStates(sampleLogs, chainStates=chainStates0)
        expect_that( length(chainStates), equals(getNChain(sDim)))
        iChain <- 1L
        chainState <- chainStates[[iChain]]
        expect_that( getChainStateParameters(chainState), equals(fx$x0Chains[,1]*iChain*iSample))
        expect_that( getLogDensityComponents(chainState), equals(logDensityComponents0*iChain*iSample))
        expect_that( isChangedByLastUpdate(chainState,blockIndices), equals(rep(FALSE,nrow(fx$x0Chains))) )
        #expect_true( is.null(getBlockIntermediate(chainState, blockIndices)) )
        iChain <- 2L    # 3L is modulo nPop equal to 1
        chainState <- chainStates[[iChain]]
        expect_that( getChainStateParameters(chainState), equals(fx$x0Chains[,1]*iChain*iSample))
        expect_that( getLogDensityComponents(chainState), equals(logDensityComponents0*iChain*iSample))
        expect_that( isChangedByLastUpdate(chainState,blockIndices), equals(rep(FALSE,nrow(fx$x0Chains))) )
        #expect_true( is.null(getBlockIntermediate(chainState, blockIndices)) )
    })


