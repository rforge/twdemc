#require(testthat)
context("DEJumpProposer")

setClass("AcceptanceTrackerMock"
        ,contains="AcceptanceTracker",
        ,representation=representation(acceptanceRate="numeric")
        ,prototype=prototype(acceptanceRate=0.25)
)
setMethod("computePopulationAcceptanceRates", signature=c("AcceptanceTrackerMock","missing"),
        function(object
                ### Average Acceptance rates within across chains.
        ){
            rep(object@acceptanceRate, getNPopulation(object@sampleDimensions) )
        }
)


# setup a SammpleLogChain
thin <- 4L
#fx <- .fixtureChainSampler(
#    nSample = c(24L,32L)
#    ,nChainInPopulation = 4L
#    ,thin=thin
#)
#sampleLogs <- .createMockSampleLogs(fx, nSampleInitial=0L)
fx=.fixtureChainSampler()
fxs <-.fixtureSampleLogsMetr(
    fx=fx
    ,sampleDimensions=.fixtureSampleDimensions(
        blockSpecifications=fx$blocks
        ,nSample = c(24L,32L)
        ,nChainInPopulation = 4L
    )
)
sampleLogs <- fxs$sampleLogs
sDim <- getSampleDimensions(sampleLogs)

acceptanceTrackerMock <- resetAcceptanceTracker( new("AcceptanceTrackerMock")
    ,sDim
    ,nSample=getNSamplePopulations(fxs$sampleDimensions)
    ,nIntervalPerRange = 2L
)
#computePopulationAcceptanceRates(acceptanceTrackerMock)

acceptanceTrackerMockLow <- resetAcceptanceTracker( new("AcceptanceTrackerMock", acceptanceRate=0.08)
        ,sDim
        ,nSample=getNSamplePopulations(fxs$sampleDimensions)
        ,nIntervalPerRange = 2L
)





jumpProposer <- initializeDimensionalSettings( new("DEJumpProposer")
    ,sampleDimensions = sDim
    ,thin=thin
    ,iParametersThatRequireJumpProposal = 2:5   # without "a"
)


test_that(".sampleStates",{
        mZ <- 4L
        #fs <- .fixtureSampleLogs()
        parametersChains <- getParametersForPopulation(sampleLogs,1L)
        res <- .sampleStates(parametersChains , mZ=mZ, nStepBack=3L, nSample=2L)
        expect_that( dim(res), equals(c( getNParameter(sDim), getNChainInPopulation(sDim)*2L, 4)))  # y:2 samples per chain, z:4 vectors for each jump
        #parametersChains <- getParameters(fs$sampleLogs[[1]])
        expect_that( res[,c(1:4),4], equals(parametersChains[,mZ,]) ) # current state, the first step, 4 chains
        expect_that( res[,c(5:8),4], equals(parametersChains[,mZ,]) ) # current state, the second step, 4 chains
        #expect_that( res[,c(1:2),4], equals(parametersChains[,mZ,]) ) # current state, the first step, 2 chains
        #expect_that( res[,c(3:4),4], equals(parametersChains[,mZ,]) ) # current state, the second step
        expect_true( all(res[c("a","m","c"),,] == 0) )    # those components not changed by multiplication
        expect_true( all(res["isM",,] == res["isC",,]) )  # 1 multiplied by same factor, difference vectors also equal
        #expect_true( all(res[c("isM","isC"),,1] != res[c("isM","isC"),,2]) )    # vector 1 disticnt from vector 2
        #expect_true( all(res[c("isM","isC"),,3] != res[c("isM","isC"),,4]) )    # vector 3 disticnt from initial state
    })

test_that("proposeJumps",{
        nGeneration <-2L*thin 
        jumpProposer <- adjustToAcceptanceRate(jumpProposer, acceptanceTracker = acceptanceTrackerMock )
        res <- proposeJumps( jumpProposer
            , nGeneration
            , sampleLogs=sampleLogs
            , iPopulationsStep= 1:getNPopulation(sDim)
            , iCurrentSamples = getNSamplePopulations(sDim)
        )
        nParameter <- getNParameter(sampleLogs) 
        nChain <- getNChain(sampleLogs)
        expect_that( dim(res$jump) , equals(c(nParameter, nGeneration, nChain)))
        expect_true( all(res$jump[c("a","m","c"),,] == 0) )   # those components all in zero state with test setup -> no differences
        expect_true( all(res$jump[c("isM","isC"),,] != 0) )   # real differences
        iChain <-1
        jumpChain <- res$jump[,,iChain]
        #
        rExtra <- res$logSnookerDensityMultiplier
        expect_that( dim(rExtra) , equals(c(nGeneration, getNChain(sDim))))
        # MAYBE: study range of rExtra
    })

test_that("proposeJumps several groups",{
            
            nGeneration <-2L*thin 
            jumpProposer <- adjustToAcceptanceRate(jumpProposer, acceptanceTracker = acceptanceTrackerMock )
            res <- proposeJumps( jumpProposer
                    , nGeneration
                    , sampleLogs=sampleLogs
                    , iPopulationsStep= 1:getNPopulation(sDim)
                    , iCurrentSamples = getNSamplePopulations(sDim)
            )
            nParameter <- getNParameter(sampleLogs) 
            nChain <- getNChain(sampleLogs)
            expect_that( dim(res$jump) , equals(c(nParameter, nGeneration, nChain)))
            expect_true( all(res$jump[c("a","m","c"),,] == 0) )   # those components all in zero state with test setup -> no differences
            expect_true( all(res$jump[c("isM","isC"),,] != 0) )   # real differences
            iChain <-1
            jumpChain <- res$jump[,,iChain]
            #
            rExtra <- res$logSnookerDensityMultiplier
            expect_that( dim(rExtra) , equals(c(nGeneration, getNChain(sDim))))
            # MAYBE: study range of rExtra
        })


test_that("proposeJumps single populations",{
        nGeneration <-2L*thin 
        jumpProposer <- adjustToAcceptanceRate(jumpProposer, acceptanceTracker = acceptanceTrackerMock )
        res <- proposeJumps( jumpProposer
            , nGeneration
            , sampleLogs=sampleLogs
            , iPopulationsStep= 2L
            , iCurrentSamples = getNSamplePopulations(sDim)
        )
        nParameter <- getNParameter(sampleLogs) 
        nChainInPopulation <- getNChainInPopulation(sampleLogs)
        expect_that( dim(res$jump) , equals(c(nParameter, nGeneration, nChainInPopulation)))
        expect_true( all(res$jump[c("a","m","c"),,] == 0) )   # those components all in zero state with test setup -> no differences
        expect_true( all(res$jump[c("isM","isC"),,] != 0) )   # real differences
        iChain <-1
        jumpChain <- res$jump[,,iChain]
        #
        rExtra <- res$logSnookerDensityMultiplier
        expect_that( dim(rExtra) , equals(c(nGeneration, nChainInPopulation)))
        # MAYBE: study range of rExtra
    })

test_that("nGenBack maxInt if not in burnin",{
            j2 <- jumpProposer
            j2@isBurnin <- FALSE
            nPop <- getNPopulation(sampleLogs)
            expect_that( .getNStepBack(j2, rep(0.25, nPop)), equals(rep(.Machine$integer.max,nPop)) )
        })

test_that("adjustToAcceptanceRate",{
            jumpProposer2 <- adjustToAcceptanceRate(jumpProposer, acceptanceTracker = acceptanceTrackerMockLow )
            expect_that( jumpProposer2@gammaAdjustmentPops[1], equals(1*0.9) )
            jumpProposer3 <- adjustToAcceptanceRate(jumpProposer2, acceptanceTracker = acceptanceTrackerMockLow )
            expect_that( jumpProposer3@gammaAdjustmentPops[1], equals(1*0.9*0.9) )
            #
            jumpProposerNoAdjustment <- setDeSetting(jumpProposer, "isAdjustingStepLengthToAcceptanceRate", FALSE )
            expect_that( getDeSettings(jumpProposerNoAdjustment)["isAdjustingStepLengthToAcceptanceRate"], equals(FALSE) )
            jumpProposer4 <- adjustToAcceptanceRate(jumpProposerNoAdjustment, acceptanceTracker = acceptanceTrackerMockLow )
            expect_that( jumpProposer4@gammaAdjustmentPops[1], equals(1) )
        })



