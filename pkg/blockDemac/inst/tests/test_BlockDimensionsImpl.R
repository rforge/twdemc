#require(testthat)
context("BlockDimensionsImpl")

x0 <- c(a=0, b=0, c=0)

blocks <- list(
    obs1 = blockSpecMock("a",, "obs1"), 
    obs2parms = blockSpecMock(c("b","c"),,c("obs2","parms"))
)

test_that("initialization",{
        blockDimensions <- newBlockDimensions( blocks, names(x0) )
        expect_that( getParameterNames(blockDimensions), equals(c("a","b","c")) )        
        expect_that( getIParametersUpdatedByBlock(blockDimensions,2L), equals(c(b=2,c=3)) )        
        #expect_that( getIParametersUsedByBlock(blockDimensions,"obs2parms"), equals(c(a=1,b=2,c=3)) )        
        expect_that( getLogDensityComponentNames(blockDimensions), equals(c("obs1","obs2","parms")) )        
        expect_that( getILogDensityComponentsByBlock(blockDimensions,2L), equals(c(obs2=2,parms=3)) )        
        expect_that( getILogDensityComponentsByBlock(blockDimensions,"obs2parms"), equals(c(obs2=2,parms=3)) )        
    })

test_that("error on non-unique blocknames",{
        names(blocks) <- NULL
        expect_error(
            blockDimensions <- newBlockDimensions( blocks, names(x0) )
        )
        names(blocks) <- rep("obs1", length(blocks))
        expect_error(
            blockDimensions <- newBlockDimensions( blocks, names(x0) )
        )
    })

test_that("correctly filled omittings in specification",{
        blocks <- list(
            obs1 = blockSpecMock("a",, "obs1"), 
            obs2parms = blockSpecMock(,c("b","c"),c("obs2","parms"))
        )
        blockDimensions <- newBlockDimensions( blocks, names(x0) )
        expect_that( getIParametersUpdatedByBlock(blockDimensions,2L), equals(c(b=2,c=3)) )        
    })


