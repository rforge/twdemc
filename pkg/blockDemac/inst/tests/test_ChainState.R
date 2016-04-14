#require(testthat)
context("ChainState")
# Furher context: fCheckColNames, fCheckColNumeric, fCheckOutsideRange

parameters = c(a=3,b=5)    
logDenComp = c(obs=12,parameters=3)
intermediates = list(results1=1:2, results2=3:4) 


chainState <- new("ChainState", 
    parameters=c(a=1,b=3,c=7), 
    logDensityComponents=c(obs1=5, obs2=7))


blockIndices <- new("BlockIndices", 
    iParametersUsed=1:2, 
    iParametersToUpdate=2L, 
    iLogDensityComponents=2L)


test_that("initialization without intermediate",{
    chainState = new("ChainState", parameters=parameters, logDensityComponents=logDenComp)
    expect_that( chainState@parameters, equals(parameters) )
    expect_that( length(chainState@isChangedByLastUpdate), equals(length(parameters)) )
    expect_that( chainState@logDensityComponents, equals(logDenComp) )
    expect_that( length(chainState@intermediates), equals(0) )
})

test_that("initialization without logDensityComponents yield empty vector",{
        chainState = new("ChainState", parameters=parameters)
        expect_that( chainState@logDensityComponents, equals(numeric(0)) )
    })

test_that("initialization without proper parameter vector throws error",{
        expect_error(
            chainState <- new("ChainState")
        )
        expect_error(
            chainState <- new("ChainState", parameters=1:3) # parameters unnamed
        )
    })
 

test_that("initialization takes intermediate",{
        chainState = new("ChainState", parameters=parameters, logDensityComponents=logDenComp, intermediates=intermediates)
        expect_that( chainState@intermediates, equals(intermediates) )
    })

test_that("initialization with constructing isChangedByLastUpdate",{
        chainState = new("ChainState", parameters=parameters )
        expect_that( chainState@isChangedByLastUpdate, equals(rep(FALSE,2)) )
        chainState = new("ChainState", parameters=parameters, isChangedByLastUpdate=TRUE )
        expect_that( chainState@isChangedByLastUpdate, equals(rep(TRUE,2)) )
        expect_error(
            chainState <- new("ChainState", parameters=parameters, isChangedByLastUpdate=c(TRUE,TRUE,FALSE) )
        )
    })

test_that("show standard",{
        chainState = new("ChainState", parameters=parameters, logDensityComponents=logDenComp, intermediates=intermediates)
        output =capture.output(chainState)
        expected <- "ChainState: parms(2)(a=3,b=5) logDen(2)(obs=12,parameters=3) interm(2)(results1,results2)\\n"
        test_that( output, equals(expected))
    })

test_that("show longParms,Den,Intermediate is shortened",{
        chainState = new("ChainState", parameters=rep(parameters,3), logDensityComponents=rep(logDenComp,3)
            , intermediates=rep(intermediates,3) )
        output <- capture.output( show(chainState) )
        expected <- "ChainState: parms(6)(a=3,b=5,a=3,b=5) logDen(6)(obs=12,parameters=3,obs=12,parameters=3) interm(6)(results1,results2,results1)\\n"
        test_that( output, equals(expected))
    })

test_that("invalidateBlockLogDensities",{
        updatedChainState <- invalidateBlockLogDensities(chainState,blockIndices) 
        expect_true( all(is.na(getBlockLogDensityComponents(updatedChainState, blockIndices))) )
    })


#--------------- getting/setting slices for specific block
test_that("getting blocks density components",{
        expect_that( getBlockLogDensityComponents(chainState, blockIndices), equals(c(obs2=7)) )
    })

test_that("getting blocks parameters",{
        expect_that( getBlockParameters(chainState, blockIndices), equals(c(a=1,b=3)) )
    })

test_that("setting blocks density components",{
        blockLogDensityComponents(chainState,blockIndices) <- 14 
        expect_that( getBlockLogDensityComponents(chainState, blockIndices), equals(c(obs2=14)) )
    })

test_that("setting blocks parameters to update",{
            # here, the mock will not invalidate other dependents
            blockUpdatersMock <- newBlockUpdaters(list(block1=blockSpecMock(names(parameters))), list(), names(parameters))
            blockParametersToUpdate(chainState, blockIndices, blockUpdatersMock) <- c(13) 
            expect_that( getBlockParameters(chainState, blockIndices), equals(c(a=1,b=13)) )
            expect_that( getBlockLogDensityComponents(chainState, blockIndices), is_equivalent_to(c(NA_real_)) )
        })

test_that("setting and getting updatedFlags",{
        isChangedByLastUpdate(chainState, blockIndices) <- TRUE
        expect_that( isChangedByLastUpdate(chainState, blockIndices), is_equivalent_to(TRUE) )
    })

test_that("setting and getting intermediates without Id defined get NULL-object: empty list",{
        expect_warning( res <- getChainStatesIntermediate(chainState, "noSet") )
        expect_that( length(res), equals(0) )
    })

test_that("setting and getting intermediates with Id defined",{
        intermediate <- list(msg="intermediate set")
        intermediateId <- "int1"
        .chainStatesIntermediate(chainState, intermediateId) <- intermediate
        expect_that( getChainStatesIntermediate(chainState, intermediateId), equals(intermediate) )
        chainState <- .invalidateChainStatesIntermediates(chainState, intermediateId)
        expect_that( length(getChainStatesIntermediate(chainState, intermediateId)), equals(0L) )
    })


