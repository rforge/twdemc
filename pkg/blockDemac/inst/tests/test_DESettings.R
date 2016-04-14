#require(testthat)
context("DESettings")

test_that("initialization empty",{
            deSettings <- new("DESettings")
            expect_that( deSettings["pSnooker"], equals(0.1))
        })

test_that("set configurable entry",{
            deSettings <- new("DESettings")
            deSettings["pSnooker"] <- 0.2
            expect_that( deSettings["pSnooker"], equals(0.2))
        })

test_that("error on setting non-configurable entry",{
            deSettings <- new("DESettings")
            expect_error(
                deSettings["nonExisting"] <- 0.2
            )
        })

test_that("initialization specific",{
            deSettings <- new("DESettings", pSnooker=0.2)
            expect_that( deSettings["pSnooker"], equals(0.2))
        })

test_that("initializeDimensionalSettings",{
            deSettings <- new("DESettings")
            expect_error(
                    deSettings["nRepresentativeRows"] <- 0.2
            )
            deSettings <- initializeDimensionalSettings(deSettings, nParameterWithProposal=5, nChainInChainGroup=4, thin=4)
            nReprRows <-  computeNRepresentativeRows(5,4)           
            expect_that( deSettings["nRepresentativeRows"], equals(nReprRows))
        })


