#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++ Unit tests for SubSpace object 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: TW
#require(testthat)
context("SubSpace")
# Furher context: fCheckColNames, fCheckColNumeric, fCheckOutsideRange

test_that("initialization",{
    subSpace = new("SubSpace")
    expect_that( length(subSpace@lowerParBounds), equals(0) )
})

test_that("split unsplitted subspace",{
        subSpace = new("SubSpace")
        subSpaces <- getSplitSubSpaces(subSpace, c(a=10) )
        expect_that( names(subSpaces), equals(c("lower","upper")))
        lower <- subSpaces$lower
        expect_that( length(getLowerParBounds(lower)), equals(0) )
        expect_true( "a" %in% names(getUpperParBounds(lower)) )
        expect_that( getUpperParBounds(lower)["a"], is_equivalent_to(10) )
        upper <- subSpaces$upper
        expect_true( "a" %in% names(getLowerParBounds(upper)) )
        expect_that( getLowerParBounds(upper)["a"], is_equivalent_to(10) )
        expect_that( length(getUpperParBounds(upper)), equals(0) )
    })


test_that("isInSubspace",{
            subSpace <- getSplitSubSpaces(new("SubSpace"), c(a=10) )$lower
            subSpace <- getSplitSubSpaces(subSpace, c(b=8) )$lower
            expect_that( getUpperParBounds(subSpace), equals(c(a=10,b=8)) )
            expect_true( isInSubSpace(subSpace, c(a=4,b=5)) )
            expect_true( isInSubSpace(subSpace, c(a=4)) )
            expect_false( isInSubSpace(subSpace, c(a=12)) )
            expect_false( isInSubSpace(subSpace, c(a=5, b=9)) )
        })

test_that("isInSubspace lower",{
            subSpace <- getSplitSubSpaces(new("SubSpace"), c(a=10) )$upper
            subSpace <- getSplitSubSpaces(subSpace, c(b=8) )$upper
            expect_that( getLowerParBounds(subSpace), equals(c(a=10,b=8)) )
            expect_true( isInSubSpace(subSpace, c(a=14,b=15)) )
            expect_true( isInSubSpace(subSpace, c(a=14)) )
            expect_false( isInSubSpace(subSpace, c(a=8)) )
            expect_false( isInSubSpace(subSpace, c(a=15, b=7)) )
        })


test_that("reflectBoundaries lower",{
            subSpace <- getSplitSubSpaces(new("SubSpace"), c(a=10) )$upper
            subSpace <- getSplitSubSpaces(subSpace, c(b=8) )$upper
            expect_that( reflectBoundaries(subSpace, c(a=14,b=15)), equals(c(a=14,b=15)) )
            expect_that( reflectBoundaries(subSpace, c(a=14)), equals(c(a=14)) )
            expect_that( reflectBoundaries(subSpace, c(a=8)), equals(c(a=12)) )
            expect_that( reflectBoundaries(subSpace, c(a=15, b=7)), equals(c(a=15, b=9)) )
        })

