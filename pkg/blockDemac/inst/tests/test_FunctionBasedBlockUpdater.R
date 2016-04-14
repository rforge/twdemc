#require(testthat)
context("FunctionBasedBlockUpdater")

chainState <- new("ChainState", parameters=c(a=1,b=3))

blockIndices <- new("BlockIndices"
        ,iParametersToUpdate = 1L
        ,iParametersUsed = 1L
        ,iParametersToUpdateInBlock = 1L
)

stepInfo <- new("StepInfoImpl", step=c(a=3,b=5))

fUpdateBlockAddTerm <- function(
        x, iParametersToUpdate
        , lowerParBounds, upperParBounds
        , intermediates
        , term
){
    list(
            isUpdated = TRUE 
            ,xC = x[ iParametersToUpdate ] + term
    )
}


test_that("initialization empty",{
            # TODO should throw an error, but validity function not called
            # validity function is never called
            generalBlockUpdater <- new("FunctionBasedBlockUpdater")
            generalBlockUpdater
            getSubSpace(generalBlockUpdater)
        })

test_that("initialization with fUpdateBlock",{
            argsF <- list(a=3)
            generalBlockUpdater <- new("FunctionBasedBlockUpdater", fUpdateBlock=identity, argsFUpdateBlock=argsF)
            expect_true( !is.null(generalBlockUpdater@fUpdateBlock) )
            expect_that( generalBlockUpdater@argsFUpdateBlock, equals(argsF) )
        })

test_that("updateParameterBlock to new chainState",{
            argsF <- list(term=2)
            generalBlockUpdater <- new("FunctionBasedBlockUpdater", fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=argsF)
            blockIndices(generalBlockUpdater) <- blockIndices
            subSpace(generalBlockUpdater) <- initializeIndexedBoundsForParameterNames(new("SubSpace"), names(chainState@parameters))
            # here, the mock will not invalidate other dependents
            blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(chainState@parameters))), list(), names(chainState@parameters))
            newChainState <- updateBlockInChainState( chainState, generalBlockUpdater, blockUpdatersMock, stepInfo )
            expect_that( newChainState@parameters, equals(c(a=3,b=3)) )
            expect_that( isChangedByLastUpdate(newChainState, generalBlockUpdater), equals(TRUE) )
        })

test_that("updateParameterBlock outside subspace range",{
            argsF <- list(term=2)
            generalBlockUpdater <- new("FunctionBasedBlockUpdater", fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=argsF)
            subSpace(generalBlockUpdater) <- getSplitSubSpaces(new("SubSpace"),c(a=2))$lower
            subSpace(generalBlockUpdater) <- initializeIndexedBoundsForParameterNames(getSubSpace(generalBlockUpdater), names(chainState@parameters))
            blockIndices(generalBlockUpdater) <- blockIndices
            blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(chainState@parameters))), list(), names(chainState@parameters))
            # note that need to assign result again to object
            expect_error(
                # update a from 1 to 3 transcends subspace boundaries
                 newChainState <- updateBlockInChainState( chainState, generalBlockUpdater, blockUpdatersMock, stepInfo )
            )
            #expect_that( isChangedByLastUpdate(newChainState, generalBlockUpdater), equals(FALSE) )
            #expect_that( newChainState@parameters, equals(c(a=1,b=3)) )
        })


#test_that("setting ChainState also changes isChanged flag",{
#            blockUpdater <- new("BlockUpdater", isChangedSinceSettingChainState = TRUE)
#            expect_that( isChangedSinceSettingChainState(blockUpdater), equals(TRUE) )
#            chainState(blockUpdater) <- chainState
#            expect_that( getChainState(blockUpdater)@parameters, equals(chainState@parameters) )
#            expect_that( isChangedSinceSettingChainState(blockUpdater), equals(FALSE) )
#        })
    