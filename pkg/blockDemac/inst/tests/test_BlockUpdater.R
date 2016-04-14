#require(testthat)
context("BlockUpdater_and_BlockIndices")

test_that("initialization empty",{
        blockIndices <- new("BlockIndices")
        blockIndices
        expect_that( length(getBlockIndicesIParametersToUpdate(blockIndices)), equals(0) )
        expect_that( length(getBlockIndicesIParametersUsed(blockIndices)), equals(0) )
        expect_that( length(getBlockIndicesILogDensityComponents(blockIndices)), equals(0) )
        expect_that( length(getBlockIndicesIntermediateIdsUsed(blockIndices)), equals(0) )
    })



context("BlockUpdater")

chainState <- new("ChainState", 
    parameters=c(a=1,b=3,c=7), 
    logDensityComponents=c(obs1=5, obs2=7))


stepInfo <- new("StepInfoImpl", step=c(a=2,b=6,c=14), temperature=c(2,3,7))

test_that("initialization empty",{
        blockUpdater <- new("BlockUpdater")
        blockUpdater
        expect_true( !is.null(getSubSpace(blockUpdater)) )
    })

test_that("initialization with subSpace",{
        splitPoint <- c(a=5)
        subSpace <- getSplitSubSpaces(new("SubSpace"),splitPoint)$lower
        blockUpdater <- new("BlockUpdater", subSpace=subSpace)
        subSpace(blockUpdater) <- subSpace
        blockUpdater
        expect_that( getUpperParBounds(getSubSpace(blockUpdater)), equals(splitPoint) )
    })




