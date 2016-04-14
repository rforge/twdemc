#require(testthat)
context("MetropolisBlockUpdater")

chainState <- new("ChainState", 
        parameters=c(a=1,b=3,c=7), 
        logDensityComponents=c(obs1=-3, obs2=-10, parms=-1)
)

blockIndices <- new("BlockIndices"
        ,iParametersToUpdate = 2L
        ,iParametersUsed = 1:2
        ,iLogDensityComponents = 2L
)

fLogDenConst <- function(
        x
        ,intermediates=list()
        ,logDensityComponents=c(obs2=NA)       ##<< numeric vector (nLDensityComponents): 
        ## already known logDensitys. Only calculate the positions that are NA
        ,...
        ,returnComponents=c(obs2=-10)
        ,logEnv=NULL        # usually in caller: new.env(parent = emptyenv())
){
    # append current inputs to log
    # its an environment, so that it can be inspected from the test outside
    if( is.environment(logEnv) ){
        logEnv$log <- paste0(logEnv$log,",x=",catNamedVector(x,3)," logDenC=",catNamedVector(logDensityComponents))
    }
    if( !is.na(logDensityComponents) ) return(logDensityComponents )
    return(returnComponents)
}


test_that("initialization with fLogDen",{
            blockUpdater <- new("MetropolisBlockUpdater"
                    , fLogDensity=fLogDenConst
                    , logDensityComponentNames = c("obs","parms")
            )
            expect_that( getArgsFLogDensity(blockUpdater), equals(list()) )
            expect_that( getMaxLogDensity(blockUpdater), equals(rep(Inf,2)) )
            #expect_true( !is.null(metropolisUpdater@fLogDensity) )
            #expect_that( metropolisUpdater@argsFLogDensity, equals(argsF) )
        })

test_that("initialization failing with not specifying logDensityComponentNames",{
            # sadly have to allow it for creating subclasses
            #expect_error(
            blockUpdater <- new("MetropolisBlockUpdater")
            #)
            expect_error(
                    blockUpdater <- new("MetropolisBlockUpdater",fLogDen=fLogDenConst)
            )
            expect_error(
                    blockUpdater <- new("MetropolisBlockUpdater",logDensityComponentNames=c("a"))
            )
        })

test_that("computeChainStatesLogDensityComponents",{
            logEnv = new.env(parent = emptyenv())
            blockUpdater <- new("MetropolisBlockUpdater"
                    , fLogDensity=fLogDenConst
                    , argsFLogDensity=list( returnComponents=c(obs2=-12), logEnv=logEnv )
                    , logDensityComponentNames = c("obs2")
                    #, blockIndices = blockIndices
            )
            blockIndices(blockUpdater) <- blockIndices
            # test logDensity not updated if chainState is not NA
            updatedChainState <- computeChainStatesLogDensityComponents(chainState, blockUpdater)
            logDenCompCurrent <- getBlockLogDensityComponents(chainState, blockUpdater)
            logDenCompNew <- getBlockLogDensityComponents(updatedChainState, blockUpdater)
            expect_that( logDenCompNew, equals(logDenCompCurrent) )
            # when chainstate parameters are set, then logDensity is set to NA, then logDensity is calculated
            # here, the mock will not invalidate other dependents
            blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(chainState@parameters))), list(), names(chainState@parameters))
            blockParametersToUpdate(updatedChainState, blockUpdater, blockUpdatersMock) <- getChainStateParameters(chainState)[ getBlockIndicesIParametersToUpdate(blockUpdater)]
            expect_true( all( is.na(getBlockLogDensityComponents(updatedChainState, blockUpdater)))) 
            updatedChainState <- computeChainStatesLogDensityComponents(updatedChainState, blockUpdater)
            logDenCompNew <- getBlockLogDensityComponents(updatedChainState, blockUpdater)
            expect_that( logDenCompNew, equals(c(obs2=-12)) )
        })


.tmp.notImplementedBecauseOfPerformance <- function(){
    test_that(".metropolisUpdatersGetProposedParametersInChainState stops on wrong step dimension",{
                blockUpdater <- new("MetropolisBlockUpdater"
                        , fLogDensity=fLogDenConst
                        , logDensityComponentNames = names(fLogDenConst(0))
                #, blockIndices = blockIndices
                )
                blockIndices(blockUpdater) <- blockIndices
                # note that need to assign result again to object
                stepInfo <- new("StepInfoImpl", step = c(a=2, b=5) )    # c is missing
                # here, the mock will not invalidate other dependents
                blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(chainState@parameters))), list(), names(chainState@parameters))
                expect_error(
                        proposedChainState <- .metropolisUpdatersGetProposedParametersInChainState(blockUpdater, chainState,  blockUpdatersMock , getStep(stepInfo) )
                )
            })
}

test_that(".metropolisUpdatersGetProposedParametersInChainState returns only parameters to update",{
            blockUpdater <- new("MetropolisBlockUpdater"
                    , fLogDensity=fLogDenConst
                    , logDensityComponentNames = names(fLogDenConst(0))
            #, blockIndices = blockIndices
            )
            blockIndices(blockUpdater) <- blockIndices
            subSpace(blockUpdater) <- initializeIndexedBoundsForParameterNames(getSubSpace(blockUpdater), names(chainState@parameters))
            # note that need to assign result again to object
            stepInfo <- new("StepInfoImpl", step = c(a=2, b=5, c=9) )    # c is missing
            # here, the mock will not invalidate other dependents
            blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(chainState@parameters))), list(), names(chainState@parameters))
            proposedChainState <- .metropolisUpdatersGetProposedParametersInChainState(blockUpdater, chainState,  blockUpdatersMock, getStep(stepInfo) )
            expect_that( getChainStateParameters(proposedChainState), equals(c(a=1,b=8,c=7)) )
        })


test_that("updateParameterBlock to new chainState",{
            logEnv = new.env(parent = emptyenv())
            blockUpdater <- new("MetropolisBlockUpdater"
                    , fLogDensity=fLogDenConst
                    , argsFLogDensity=list( returnComponents=c(obs2=-8), logEnv=logEnv )
                    , logDensityComponentNames = names(fLogDenConst(0))
            #, blockIndices = blockIndices
            )
            blockIndices(blockUpdater) <- blockIndices
            subSpace(blockUpdater) <- initializeIndexedBoundsForParameterNames(getSubSpace(blockUpdater), names(chainState@parameters))
            stepInfo <- new("StepInfoImpl", step = c(a=2, b=5, c=9), nLogDensityComponents=3 )
            # here, the mock will not invalidate other dependents
            blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(chainState@parameters))), list(), names(chainState@parameters))
            updatedChainState <- updateBlockInChainState( chainState, blockUpdater, blockUpdatersMock, stepInfo )
            expect_that( getChainStateParameters(updatedChainState), equals(c(a=1,b=8,c=7)) )
            expect_that( isChangedByLastUpdate(updatedChainState, blockUpdater), is_equivalent_to(TRUE) )
        })

test_that("updateParameterBlock outside subspace range reflected and accepted",{
            logEnv = new.env(parent = emptyenv())
            blockUpdater <- new("MetropolisBlockUpdater"
                    , fLogDensity=fLogDenConst
                    , argsFLogDensity=list( returnComponents=c(obs2=-8), logEnv=logEnv )
                    , logDensityComponentNames = names(fLogDenConst(0))
            #, blockIndices = blockIndices
            )
            blockIndices(blockUpdater) <- blockIndices
            subSpace(blockUpdater) <- getSplitSubSpaces(new("SubSpace"),c(b=7))$lower
            subSpace(blockUpdater) <- initializeIndexedBoundsForParameterNames(getSubSpace(blockUpdater), names(chainState@parameters))
            #
            stepInfo <- new("StepInfoImpl", step = c(a=2, b=5, c=9), nLogDensityComponents=3 )
            # here, the mock will not invalidate other dependents
            blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(chainState@parameters))), list(), names(chainState@parameters))
            updatedChainState <- updateBlockInChainState( chainState, blockUpdater, blockUpdatersMock, stepInfo )
            expect_that( updatedChainState@parameters, equals(c(a=1,b=6,c=7)) )
            expect_that( isChangedByLastUpdate(updatedChainState, blockUpdater), is_equivalent_to(TRUE) )
        })

test_that("updateParameterBlock not accepted",{
            logEnv = new.env(parent = emptyenv())
            blockUpdater <- new("MetropolisBlockUpdater"
                    , fLogDensity=fLogDenConst
                    , argsFLogDensity=list( returnComponents=c(obs2=-.Machine$integer.max), logEnv=logEnv )
                    , logDensityComponentNames = names(fLogDenConst(0))
            #, blockIndices = blockIndices
            )
            blockIndices(blockUpdater) <- blockIndices
            subSpace(blockUpdater) <- initializeIndexedBoundsForParameterNames(getSubSpace(blockUpdater), names(chainState@parameters))
            stepInfo <- new("StepInfoImpl", step = c(a=2, b=5, c=9), nLogDensityComponents=3 )
            # here, the mock will not invalidate other dependents
            blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(chainState@parameters))), list(), names(chainState@parameters))
            updatedChainState <- updateBlockInChainState( chainState, blockUpdater, blockUpdatersMock,stepInfo )
            expect_that( updatedChainState, equals(chainState) )
            expect_that( isChangedByLastUpdate(updatedChainState, blockUpdater), is_equivalent_to(FALSE) )
        })


