#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++ Unit tests for StepInfoRangeImpl data object 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: TW
#require(testthat)
context("StepInfoRangeImpl")
# Furher context: fCheckColNames, fCheckColNumeric, fCheckOutsideRange


blockIndices <- new("BlockIndices"
    ,iParametersToUpdate = 2L
    ,iParametersUsed = 1:2
    ,iLogDensityComponents = 2L
)
blockIndices1 <- new("BlockIndices"
        ,iLogDensityComponents = 1L
)
blockIndices2 <- new("BlockIndices"
        ,iLogDensityComponents = 1:2
)
blockIndices3 <- new("BlockIndices"
        ,iLogDensityComponents = 1:3
)


test_that("empty constructor allowed for subclasses",{
        stepInfo <- new("StepInfoRangeImpl")
})

test_that("missing step in constructor throws error",{
        expect_error(      
            stepInfo <- new("StepInfoRangeImpl", temperature=1)
        )
        expect_error(    
                # missing names
                stepInfo <- new("StepInfoRangeImpl", step=c(1))
        )
    })

test_that("rExtra and temperature with suitable defaults",{
        stepInfo = new("StepInfoRangeImpl", step=c(a=1))
        stepInfo
        #
        expect_true( implementsStepInfo(stepInfo) )
        #    
        expect_that(getRExtra(stepInfo), equals(0) )
        expect_that(getBlockTemperature(stepInfo,blockIndices1), equals(c(1)) )
        #
        stepInfo = new("StepInfoRangeImpl", step=c(a=1), temperature=5)
        expect_that(getRExtra(stepInfo), equals(0) )
        expect_that(getBlockTemperature(stepInfo,blockIndices1), equals(5) )
        #
        stepInfo = new("StepInfoRangeImpl", step=c(a=1), temperature=5, nLogDensityComponents=3)
        expect_that(getRExtra(stepInfo), equals(0) )
        expect_that(getBlockTemperature(stepInfo,blockIndices3), equals(rep(5,3)) )
    })

test_that("getting blocks temperature from stepInfo",{
        stepInfo = new("StepInfoRangeImpl", step=c(a=1,b=2), temperature=c(2,3) )
        expect_that( getBlockTemperature(stepInfo, blockIndices1), is_equivalent_to(c(2)) )
        expect_that( getBlockTemperature(stepInfo, blockIndices2), is_equivalent_to(c(2,3)) )
        expect_error(
         getBlockTemperature(stepInfo, blockIndices3)
        )
    })

test_that("setting blocks temperature in step info",{
        stepInfo = new("StepInfoRangeImpl", step=c(a=1,b=2), temperature=c(2,3), nGeneration=4 )
        expect_that( getBlockTemperature(stepInfo, blockIndices2), is_equivalent_to(c(2,3)) )
        .blockTemperature(stepInfo,blockIndices2) <- c(6,7) 
        expect_that( getBlockTemperature(stepInfo, blockIndices2), is_equivalent_to(c(6,7)) )
    })

test_that("getting temperature for second step",{
            stepInfo = new("StepInfoRangeImpl", step=matrix(1:12, nrow=2, dimnames=list(c("a","b"),NULL)), temperature=1:2 )
            expect_that( getGenerationIndex(stepInfo), equals(1L) )
            expect_that( getStep(stepInfo), equals(c(a=1,b=2)) )
            generationIndex(stepInfo) <- 5L
            expect_that( getGenerationIndex(stepInfo), equals(5L) )
            expect_that( getStep(stepInfo), equals(c(a=9,b=10)) )
        })


