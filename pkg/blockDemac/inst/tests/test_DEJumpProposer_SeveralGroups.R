#require(testthat)
context("DEJumpProposer_SeveralGroups")

setClass("AcceptanceTrackerMock"
        ,contains="AcceptanceTracker",
)
setMethod("computePopulationAcceptanceRates", signature=c("AcceptanceTrackerMock","missing"),
        function(object
                ### Average Acceptance rates within across chains.
        ){
            rep(0.25, getNPopulation(object@sampleDimensions) )
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
        ,nChainInPopulation = 6L
        ,nChainGroupsInPopulation = 2L
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




jumpProposer <- initializeDimensionalSettings( new("DEJumpProposer")
    ,sampleDimensions = sDim
    ,thin=thin
    ,iParametersThatRequireJumpProposal = 2:5   # without "a"
)

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




