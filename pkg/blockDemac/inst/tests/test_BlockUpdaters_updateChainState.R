#require(testthat)
context("BlockUpdaters:updateChainState")


x0 <- c(isA=1, A=3, isB=1, B=5, isC=1, C=7)

chainState <- new("ChainState", 
    parameters=x0, 
    logDensityComponents=c(obs1=-100, obs2=-200, obs3=-300, parms=-20)
)

stepInfo <- new("StepInfoImpl", step = c(isA=0, A=3, isB=0, B=5, isC=0, C=7), nLogDensityComponents=4 ) 

fLogDenConst <- function(
    x
    ,logDensityComponents       ##<< numeric vector (nLDensityComponents): 
    ## already known logDensitys. Only calculate the positions that are NA
    ,...
    ,acceptedLogDensityComponents
    ,logEnv=NULL        # usually in caller: new.env(parent = emptyenv())
){
    # if first parameter is 1 then accept, else not
    if( is.environment(logEnv) ){
        logEnv$log <- paste0(logEnv$log,",x=",catNamedVector(x,3)," logDenC=",catNamedVector(logDensityComponents))
    }
    ret <- if (x[1]==1) {
            acceptedLogDensityComponents    # accepted  
        } else {
            Inf*acceptedLogDensityComponents    # not accepted
        }
    ret
}

# due to 0 in stepInfo isA and isB stay at their original values 
blocks <- list(
    lden1 = blockSpec(,c("isA","A"), new("MetropolisBlockUpdater", 
            fLogDensity=fLogDenConst,
            argsFLogDensity=list(acceptedLogDensityComponents = c(obs1=-10)),
            logDensityComponentNames = c("obs1"))), 
    lden2 = blockSpec(c("isC","C"),names(x0),new("MetropolisBlockUpdater",    # depends on all parameters 
            fLogDensity=fLogDenConst,
            argsFLogDensity=list(acceptedLogDensityComponents = c(obs2=-20)),
            logDensityComponentNames = c("obs2"))),
    lden34 = blockSpec(,c("isB","B"),new("MetropolisBlockUpdater", 
            fLogDensity=fLogDenConst,
            argsFLogDensity=list(acceptedLogDensityComponents = c(parms=-2, obs3=-30)),
            logDensityComponentNames = c("obs3","parms")))
)

intermediateSpecs <- list()


#blockDimensions <- newBlockDimensions(blocks, names(x0))
blockUpdaters <- blockUpdatersInit <- newBlockUpdaters(blocks, intermediateSpecs, names(x0))
lden1 <- getBlockUpdaterForName(blockUpdatersInit,"lden1") 
lden2 <- getBlockUpdaterForName(blockUpdatersInit,"lden2") 
lden34 <- getBlockUpdaterForName(blockUpdatersInit,"lden34") 

test_that(".initializeBlockDependencies",{
        #blockUpdaters <- .initializeBlockIndices(blockUpdaters, names(x0) )
        #blockUpdaters <- .initializeBlockDependencies(blockUpdaters)
        lden2 <- getBlockUpdaterForName(blockUpdaters,"lden2") 
        # blocks 1 and 3 have dependent block 2,
        expect_that( length(getDependentBlockUpdaters(blockUpdaters,1L)), equals(2L) )
        expect_that( getDependentBlockUpdaters(blockUpdaters,1L)[[2]], equals(lden2) )
        expect_that( getDependentBlockUpdaters(blockUpdaters,"lden1")[[2]], equals(lden2) )
        expect_that( length(getDependentBlockUpdaters(blockUpdaters,3L)), equals(2L) )
        expect_that( getDependentBlockUpdaters(blockUpdaters,3L)[[1]], equals(lden2) )
        # block 2 has no dependent block
        expect_that( length(getDependentBlockUpdaters(blockUpdaters,2L)), equals(1L) )
    })


test_that(".invalidateBlocksThatDependOn",{
        blockDim <- getBlockDimensions(blockUpdatersInit)
        expect_that( getLogDensityComponentNames(blockDim)
            , is_equivalent_to(names(chainState@logDensityComponents)))     
        updatedChainState <- .invalidateBlocksThatDependOn(chainState, blockUpdatersInit, 1L)
        expect_true( all( is.na(getBlockLogDensityComponents(updatedChainState,lden1))) )
        expect_true( all( is.na(getBlockLogDensityComponents(updatedChainState,lden2))) )
        expect_true( all( !is.na(getBlockLogDensityComponents(updatedChainState,lden34))) )
    })


test_that("computeUpdatedChainState all accepted",{
        #stepInfo@step
        updatedChainState <- computeUpdatedChainState(chainState, blockUpdatersInit, stepInfo)
        expect_that( isChangedByLastUpdate(updatedChainState,lden1), equals(c(TRUE,TRUE)) )
        expect_that( isChangedByLastUpdate(updatedChainState,lden2), equals(c(TRUE,TRUE)) )
        expect_that( isChangedByLastUpdate(updatedChainState,lden34), equals(c(TRUE,TRUE)) )
        expect_that( getChainStateParameters(updatedChainState), equals(c(isA=1, A=6, isB=1, B=10, isC=1, C=14)))
        # lden4 invalidated by lden23        
        expect_that( getBlockLogDensityComponents(updatedChainState,lden2), equals(c(obs2=NA_real_)))
    })

test_that("computeUpdatedChainState second not accepted",{
        chainStateIsB0 <- chainState
        chainStateIsB0@parameters["isB"] <- 0   # B (lDen34) will not be updated
        updatedChainState <- computeUpdatedChainState(chainStateIsB0, blockUpdatersInit, stepInfo)
        expect_that( isChangedByLastUpdate(updatedChainState,lden1), equals(c(TRUE,TRUE)) )
        expect_that( isChangedByLastUpdate(updatedChainState,lden34), equals(c(FALSE,FALSE)) )
        expect_that( getChainStateParameters(updatedChainState), equals(c(isA=1, A=6, isB=0, B=5, isC=1, C=14)))
        # accepted: logDen from fLogDen
        expect_that( getBlockLogDensityComponents(updatedChainState,lden1), equals(c(obs1=-10)))
        # this time not invalidated by update of lden34 block (executed after lden2), because lden34 was not accepted
        expect_that( getBlockLogDensityComponents(updatedChainState,lden2), equals(c(obs2=-20)))
        # not accepted: chainState as previous 
        expect_that( getBlockLogDensityComponents(updatedChainState,lden34), equals(c(obs3=-300, parms=-20)))
    })

test_that("blockUpdaterApply",{
            chainStateIsB0 <- chainState
            chainStateIsB0@parameters["isB"] <- 0   # B will not be updated
            updatedChainState <- computeUpdatedChainState(chainStateIsB0, blockUpdatersInit, stepInfo)
            res <- blockUpdaterApply(blockUpdatersInit, function(blockUpdater, chainState=updatedChainState){
                        any(isChangedByLastUpdate( chainState, blockUpdater))
                    })   
            expect_that( res, is_equivalent_to(c(TRUE, TRUE, FALSE)))
        })







