#require(testthat)
context("InternalMetropolisBlockUpdater")

chainState <- new("ChainState", 
    parameters=c(a=1,b=3,c=7), 
    logDensityComponents=c(obs1=-3, obs2=-10, parms=-1)
)

blockIndices <- new("BlockIndices"
    ,iParametersToUpdate = 2L   # a,b
    ,iParametersUsed = 1:2      # a,b
    ,iLogDensityComponents = 2:3   # obs2, parms   # note here two 
)

fLogDenInternalConst <- function(
    # will do internal rejection if x["a"]==0 and return logDensity of obs2 as NA
    x
    ,intermediates=list()
    ,logDensityComponents=c(obs2=NA,parms=NA)	##<< scalar: logDen for parms for current x (to 
    ,logDensityComponentsAccepted=numeric(0)    ##<< scalar: logDen for parms from revious run for two step Metropolis decision
    ,temperature=1		                        ##<< numeric named vector: the temperature for internal metropolis step
    ,isInternalRejectionUsed=TRUE               ##<< set to FALSE to enforce calculation of all components
    ,...
    ,returnComponents=c(obs2=-10,parms=-1)
    ,logEnv=NULL        # usually in caller: new.env(parent = emptyenv())
){
    # append current inputs to log
    # its an environment, so that it can be inspected from the test outside
    if( is.environment(logEnv) ){
        logEnv$log <- paste0(logEnv$log,",x=",catNamedVector(x,3)," logDenC=",catNamedVector(logDensityComponents))
    }
    assert_that( !isInternalRejectionUsed || length(logDensityComponentsAccepted) )
    newLogDen = c(obs2=NA,parms=NA)
    # calculate prior (here just assign)
    newLogDen["parms"] <- returnComponents["parms"]
    # perform internal rejection
    if( isInternalRejectionUsed ){
        if( x["a"] == 0 ) return( newLogDen) 
    } 
    # calculate expensive model if not internally rejected
    newLogDen["obs2"] <- returnComponents["obs2"]
    newLogDen
}

test_that("initialization failing with not specifying logDensityComponentNames",{
        # sadly have to allow it for creating subclasses
        #expect_error(
        blockUpdater <- new("InternalMetropolisBlockUpdater")
        #)
        expect_error(
            blockUpdater <- new("InternalMetropolisBlockUpdater",fLogDen=fLogDenConst)
        )
        expect_error(
            blockUpdater <- new("InternalMetropolisBlockUpdater",logDensityComponentNames=c("a"))
        )
    })


test_that("computeChainStatesLogDensityComponents internally rejected",{
        logEnv = new.env(parent = emptyenv())
        blockUpdater <- new("InternalMetropolisBlockUpdater"
            , fLogDensity=fLogDenInternalConst
            , argsFLogDensity=list( returnComponents=c(obs2=-6,parms=-2), logEnv=logEnv )
            , logDensityComponentNames = c("obs2","parms")
    )
    blockIndices(blockUpdater) <- blockIndices
    #getChainStateParameters(chainState)
        chainStateParameters(chainState) <- c(0,3,7)
        # default without internal rejection
        newChainState <- computeChainStatesLogDensityComponents(chainState, blockUpdater)
        expect_that( getBlockLogDensityComponents(newChainState, blockUpdater), equals(c(obs2=-6,parms=-2)))
        # internal rejection
        newChainState <- computeChainStatesLogDensityComponents(chainState, blockUpdater, isInternalRejectionUsed=TRUE)
        expect_that( getBlockLogDensityComponents(newChainState, blockUpdater), equals(c(obs2=NA_real_,parms=-2)))
        #logEnv$log
    })

test_that("updateParameterBlock to new chainState",{
        logEnv = new.env(parent = emptyenv())
        blockUpdater <- new("InternalMetropolisBlockUpdater"
            , fLogDensity=fLogDenInternalConst
            , argsFLogDensity=list( returnComponents=c(obs2=-6,parms=-2), logEnv=logEnv )
            , logDensityComponentNames = c("obs2","parms")
    )
    blockIndices(blockUpdater) <- blockIndices
    subSpace(blockUpdater) <- initializeIndexedBoundsForParameterNames(getSubSpace(blockUpdater), names(chainState@parameters))
    stepInfo <- new("StepInfoImpl", step = c(a=2, b=5, c=9), temperature=2, nLogDensityComponents=3 )
        # here, the mock will not invalidate other dependents
            blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(chainState@parameters))), list(), names(chainState@parameters))
            newChainState <- updateBlockInChainState(chainState, blockUpdater,   blockUpdatersMock, stepInfo )
        expect_that( isChangedByLastUpdate(newChainState, blockUpdater), is_equivalent_to(TRUE) )
        expect_that( newChainState@parameters, equals(c(a=1,b=8,c=7)) )
        expect_that( getBlockLogDensityComponents(newChainState, blockUpdater), equals(c(obs2=-6,parms=-2)) )
    })

test_that("updateParameterBlock rejected by internal metropolis",{
        logEnv = new.env(parent = emptyenv())
        logDenCompOrig <- getLogDensityComponents(chainState)
        blockUpdater <- new("InternalMetropolisBlockUpdater"
            , fLogDensity=fLogDenInternalConst
            , argsFLogDensity=list( returnComponents=c(obs2=-6,parms=-2), logEnv=logEnv )
            , logDensityComponentNames = c("obs2","parms")
    )
    blockIndices(blockUpdater) <- blockIndices
    subSpace(blockUpdater) <- initializeIndexedBoundsForParameterNames(getSubSpace(blockUpdater), names(chainState@parameters))
    chainStateParameters(chainState) <- c(0,3,7)    
        # setting parameters invalidates logDensity, we want to check that is not replaced by fLogDen
        chainState@logDensityComponents  <- logDenCompOrig
        #blockParameters(chainState, blockUpdater) <- c(0,3)
        stepInfo <- new("StepInfoImpl", step = c(a=2, b=5, c=9), temperature=2, nLogDensityComponents=3 )
        blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(chainState@parameters))), list(), names(chainState@parameters))
        newChainState <- updateBlockInChainState( chainState, blockUpdater, blockUpdatersMock ,stepInfo )
        expect_that( isChangedByLastUpdate(newChainState, blockUpdater), is_equivalent_to(FALSE) )
        expect_that( getChainStateParameters(newChainState), equals( getChainStateParameters(chainState)) )
        expect_that( getLogDensityComponents(newChainState), equals(logDenCompOrig) )
    })





