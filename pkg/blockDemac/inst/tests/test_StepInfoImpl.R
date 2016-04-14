#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++ Unit tests for StepInfoImpl data object 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: TW
#require(testthat)
context("StepInfoImpl")
# Furher context: fCheckColNames, fCheckColNumeric, fCheckOutsideRange


blockIndices <- new("BlockIndices"
    ,iParametersToUpdate = 2L
    ,iParametersUsed = 1:2
    ,iLogDensityComponents = 2L
)


test_that("empty constructor allowed for subclasses",{
        stepInfo <- new("StepInfoImpl")
        expect_true( implementsStepInfo(stepInfo) )
        
})

test_that("missing step in constructor throws error",{
        expect_error(      
            stepInfo <- new("StepInfoImpl", temperature=1)
        )
    })

test_that("rExtra and temperature with suitable defaults",{
        stepInfo = new("StepInfoImpl", step=c(a=1))
        stepInfo
        expect_that(getRExtra(stepInfo), equals(0) )
        expect_that(stepInfo@temperature, equals(1) )
        #
        stepInfo = new("StepInfoImpl", step=c(a=1), temperature=5)
        expect_that(getRExtra(stepInfo), equals(0) )
        expect_that(stepInfo@temperature, equals(5) )
        #
        stepInfo = new("StepInfoImpl", step=c(a=1), temperature=5, nLogDensityComponents=3)
        expect_that(getRExtra(stepInfo), equals(0) )
        expect_that(stepInfo@temperature, equals(rep(5,3)) )
        
    })

test_that("getting blocks temperature from stepInfo",{
        stepInfo = new("StepInfoImpl", step=c(a=1), temperature=c(1,3,7))
        expect_that( getBlockTemperature(stepInfo, blockIndices), is_equivalent_to(3) )
    })

test_that("setting blocks temperature in step info",{
        stepInfo = new("StepInfoImpl", step=c(a=1), temperature=c(1,3,7))
        .blockTemperature(stepInfo,blockIndices) <- 6 
        expect_that( getBlockTemperature(stepInfo, blockIndices), is_equivalent_to(6) )
    })

