#require(testthat)
context("SampleDimensionsImpl")


thin <- 4L
nSamplePopulations <- c(16L,32L)
nChainInPopulation <- 4L
nPopulation <- length(nSamplePopulations)
fx <- .fixtureChainSampler(thin=thin)

sampleDimensions <- initializeSampleDimensionsImpl( new("SampleDimensionsImpl")
        , blockDimensions=fx$blockDimensions
        , nSamplePopulations=nSamplePopulations, nChainInPopulation=nChainInPopulation)

test_that("empty initialization",{
        emptyObject <- new("SampleDimensionsImpl")
    })

#test_that("initialization without blockupdaters",{
#            sampleDimensions <- new("SampleDimensionsImpl", blockUpdaters=fx$blockUpdaters, thin=thin, nSamplePopulations=nSamplePopulations, nChainInPopulation=nChainInPopulation)
#            expect_that( getNPopulation(emptyObject), equals(0L))
#        })

test_that("getSampleDimensions blocks",{
        #expect_that( getThin(sampleDimensions), equals(thin))
        nBlock <- getNBlock(fx$blockDimensions)
        expect_that( getNBlock(sampleDimensions), equals(nBlock))
        expect_that( getBlockNames(sampleDimensions), equals(c("add7", "met1", "met2")))
    })


test_that("getSampleDimensions parameters",{
        expect_that( getNParameter(sampleDimensions), equals(length(fx$x0)))
        expect_that( getParameterNames(sampleDimensions), equals(names(fx$x0)))
        expect_that( getIParametersUpdatedByBlock(sampleDimensions,2L), equals(c(isM=2,m=3)))
        expect_that( getIParametersUpdatedByBlock(sampleDimensions,"met1"), equals(c(isM=2,m=3)))
    })

test_that("getSampleDimensions logDensityComponents",{
        expect_that( getNLogDensityComponent(sampleDimensions), equals(length(fx$logDensityComponents0)))
        expect_that( getLogDensityComponentNames(sampleDimensions), equals(names(fx$logDensityComponents0)))
        expect_that( getILogDensityComponentsByBlock(sampleDimensions,3L), equals(c(m2=2,parms=3)))
        expect_that( getILogDensityComponentsByBlock(sampleDimensions,"met2"), equals(c(m2=2,parms=3)))
    })

test_that("getSampleDimensions populations",{
        expect_that( getNSamplePopulations(sampleDimensions), equals(nSamplePopulations))
        expect_that( getNPopulation(sampleDimensions), equals(nPopulation))
        expect_that( getNChainInPopulation(sampleDimensions), equals(nChainInPopulation))
        expect_that( getNChain(sampleDimensions), equals(nPopulation*nChainInPopulation))
        expect_that( getIChainsForPopulation(sampleDimensions,2L), equals(5:8))
        expect_error( getIChainsForPopulation(sampleDimensions,nPopulation+1))
        expect_that( getIChainsForPopulations(sampleDimensions,1:2), equals(1:8))
        expect_error( getIChainsForPopulations(sampleDimensions,1:3))
        expect_that( getIPopulationsForChains(sampleDimensions,1:13), equals(c(1,1,1,1,2,2,2,2,3,3,3,3,4)))
    })



