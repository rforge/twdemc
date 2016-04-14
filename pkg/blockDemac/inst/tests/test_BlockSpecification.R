#require(testthat)
context("BlockSpecification")


fUpdateBlockAddTerm <- function(
    x, iParametersToUpdate
    , lowerParBounds, upperParBounds, intermediate
    , term
){
    list(
        isUpdated = TRUE 
        ,xC = x[ iParametersToUpdate ] + term
        ,intermediate = NULL    # actually dont return NULL, better list(), here testing handled
    )
}
blockUpdater <- new("FunctionBasedBlockUpdater", 
    fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=7))

test_that("initialization",{
        # typical
        blockSpecification <- blockSpec("a",c("a","b"),blockUpdater )
        # all parameters used
        blockSpecification <- blockSpec("a",, blockUpdater)
        expect_that( getParametersToUpdate(blockSpecification) , equals(c("a")))
        expect_that( getParametersUsed(blockSpecification) , equals(c("a")))
        # additionally all parameters update
        blockSpecification <- blockSpec(,, blockUpdater)
        expect_that( getParametersToUpdate(blockSpecification) , equals(character(0)))
        expect_that( getParametersUsed(blockSpecification) , equals(character(0)))
        # only b parameter used and updated
        blockSpecification <- blockSpec(,"b",blockUpdater )
        expect_that( getParametersToUpdate(blockSpecification) , equals(c("b")) )
        expect_that( getParametersUsed(blockSpecification), equals(c("b")))
    })

test_that("missing BlockUpdater throws error",{
        expect_error(
            blockSpecification <- blockSpec("a",, NULL)
        )
    })


test_that("parmsToUpdate not in parmsUsed gives warning",{
        expect_warning(
            blockSpecification <- blockSpec("a","b",blockUpdater )
        )
    })

