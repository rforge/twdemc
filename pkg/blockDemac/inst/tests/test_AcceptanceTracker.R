#require(testthat)
context("AcceptanceTracker")

fx=.fixtureChainSampler()
fxs <-.fixtureSampleLogsMetr(
    fx=fx
    ,sampleDimensions=.fixtureSampleDimensions(
        blockSpecifications=fx$blocks
        ,nSample = rep(c(24L),2L)
        ,nChainInPopulation = 2L
    )
) 

acceptanceTracker <- new("AcceptanceTracker", settings=new("AcceptanceTrackerSettings"))
nBlock <- getNBlock(fxs$sampleDimensions) 
nSample <- max( getNSamplePopulations(fxs$sampleDimensions))
nIntervalPerRange <- 2
#iPopulationChains <- getIPopulationChains(fx$sampleDimensions)
nChain <- getNChain(fxs$sampleDimensions) #length(iPopulationChains)

test_that("resetAcceptanceTracker",{
        nSampleBatch <- nIntervalPerRange*c(5L,7L)
        acceptanceTracker <- resetAcceptanceTracker(acceptanceTracker
            ,sampleDimensions=fxs$sampleDimensions
            ,nSample = nSampleBatch
            , nIntervalPerRange=nIntervalPerRange)
        proportionAcceptedWindow <- acceptanceTracker@proportionAcceptedWindow
        proportionAcceptedWindow
        expect_true( is.array(proportionAcceptedWindow) )
        dimW <- dim(proportionAcceptedWindow)
        expect_that( dimW, is_equivalent_to(c(nBlock, 2*5+2, nChain)) )
        #
        acceptanceLog <- getAcceptanceLog(acceptanceTracker) 
        expect_true( is.array(acceptanceLog) )
        dimW <- dim(acceptanceLog)
        expect_that( dimW, is_equivalent_to(c(nBlock, max(nSampleBatch), nChain)) )
        expect_that( getCountSample(acceptanceTracker), equals(0L))        
    })

test_that("recordAcceptanceRate",{
        nSampleBatch <- nIntervalPerRange*c(5L,5L)
        acceptanceTracker <- resetAcceptanceTracker(acceptanceTracker
            ,sampleDimensions=fxs$sampleDimensions
            ,nSample = nSampleBatch
            , nIntervalPerRange=nIntervalPerRange)
        proportionAcceptedInInterval <- array(1:(nBlock*nIntervalPerRange)/2, dim=c(nBlock, nIntervalPerRange, nChain ))
        acceptanceTracker <- recordAcceptanceRate(acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval)
        alog <- getAcceptanceLog(acceptanceTracker)
        alog
        expect_true( !any(is.na(alog[,1:2,]) ))
        for( i in 1:4 ){
            acceptanceTracker <- recordAcceptanceRate(acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval)
        }
        nSample <- getCountSample(acceptanceTracker)
        expect_that( nSample, equals(2L*5L))
        alog <- getAcceptanceLog(acceptanceTracker)
        alog
        expect_true( !any(is.na(alog) ))
        expect_error(
            # not enough space in log
            acceptanceTracker <- recordAcceptanceRate(acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval)
        )
    }
)

test_that("recordAcceptanceRate_SomeChains",{
        nSampleBatch <- nIntervalPerRange*c(12L,12L)
        acceptanceTracker <- resetAcceptanceTracker(acceptanceTracker
            ,sampleDimensions=fxs$sampleDimensions
            ,nSample = nSampleBatch
            , nIntervalPerRange=nIntervalPerRange)
        iChains <- 3:4
        proportionAcceptedInInterval <- array(1:(nBlock*nIntervalPerRange)/2, dim=c(nBlock, nIntervalPerRange, length(iChains) ))
        acceptanceTracker <- recordAcceptanceRate(acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval, iChains=iChains)
        alog <- getAcceptanceLog(acceptanceTracker)
        alog
        expect_true( !any(is.na(alog[,1:2,3:4]) ))
        expect_true( all(is.na(alog[,,1:2]) ))
        for( i in 1:11 ){
            acceptanceTracker <- recordAcceptanceRate(acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval, iChains=iChains)
        }
        nSample <- getCountSample(acceptanceTracker)
        expect_that( nSample, equals(2L*12L))
        alog <- getAcceptanceLog(acceptanceTracker)
        alog
        expect_true( !any(is.na(alog[,,3:4]) ))
        expect_true( all(is.na(alog[,,1:2]) ))
        expect_error(
            # not enough space in log
            acceptanceTracker <- recordAcceptanceRate(acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval, iChains=iChains)
        )
    })


test_that("computePopulationAcceptanceRates",{
        nSampleBatch <- getNSamplePopulations(fxs$sampleDimensions)
        acceptanceTracker <- resetAcceptanceTracker(acceptanceTracker
            ,sampleDimensions=fxs$sampleDimensions
            ,nSample = nSampleBatch
            , nIntervalPerRange=nIntervalPerRange)
        popAcceptanceRates <- computePopulationAcceptanceRates(acceptanceTracker)
        expect_that( popAcceptanceRates, equals(c(0.25,0.25))) # initial acceptancerate in settings
        #
        proportionAcceptedInInterval <- array(0.5, dim=c(nBlock, nIntervalPerRange, nChain ))
        acceptanceTracker <- recordAcceptanceRate(acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval)
        popAcceptanceRates <- computePopulationAcceptanceRates(acceptanceTracker)
        pExp <- (0.25+2*0.5)/3  # two thinning intervals of 0.5 acceptance
        expect_that( popAcceptanceRates, equals(rep(pExp,2))) 
    })

test_that("computePopulationAcceptanceRates_somePopulations",{
        nSampleBatch <- getNSamplePopulations(fxs$sampleDimensions)
        acceptanceTracker <- resetAcceptanceTracker(acceptanceTracker
            ,sampleDimensions=fxs$sampleDimensions
            ,nSample = nSampleBatch
            , nIntervalPerRange=nIntervalPerRange)
        popAcceptanceRates <- computePopulationAcceptanceRates(acceptanceTracker, 2L)
        expect_that( popAcceptanceRates, equals(c(0.25))) # initial acceptancerate in settings
        #
        iChains <- 3:4
        proportionAcceptedInInterval <- array(0.5, dim=c(nBlock, nIntervalPerRange, length(iChains) ))
        acceptanceTracker <- recordAcceptanceRate(acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval, iChains=iChains)
        popAcceptanceRates <- computePopulationAcceptanceRates(acceptanceTracker,2L)
        pExp <- (0.25+2*0.5)/3
        expect_that( popAcceptanceRates, equals(rep(pExp,1))) 
    })

test_that("computeChainGroupAcceptanceRate",{
            nSampleBatch <- getNSamplePopulations(fxs$sampleDimensions)
            acceptanceTracker <- resetAcceptanceTracker(acceptanceTracker
                    ,sampleDimensions=fxs$sampleDimensions
                    ,nSample = nSampleBatch
                    , nIntervalPerRange=nIntervalPerRange)
            groupAcceptanceRates <- computeChainGroupAcceptanceRate(acceptanceTracker, 2L, 1:2)
            expect_that( groupAcceptanceRates, equals(c(0.25,0.25))) # initial acceptancerate in settings
            #
            proportionAcceptedInInterval <- array(0.5, dim=c(nBlock, nIntervalPerRange, nChain ))
            acceptanceTracker <- recordAcceptanceRate(acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval)
            groupAcceptanceRates <- computeChainGroupAcceptanceRate(acceptanceTracker, 2L, 1:2)
            pExp <- (0.25+2*0.5)/3  # two thinning intervals
            expect_that( groupAcceptanceRates, equals(rep(pExp,1L))) 
        })







