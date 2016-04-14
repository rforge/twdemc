#' @include BlockUpdater.R
#' @include BlockUpdaters.R

#' @export
setClass("FunctionBasedBlockUpdater"
    , contains ="BlockUpdater"	    # blockIndices, subSpace, chainState
    , representation(  
        fUpdateBlock="function"	    ##<< f(theta, iParms, upperParBounds, lowerParBounds, intermediate, ...) 
        , argsFUpdateBlock = "list"
    )
    ,validity=function(object){
        if( is.null(body(object@fUpdateBlock)) ) return("Need to provide an update function.")
        TRUE
    }
)


#' @export
setMethod("updateBlockInChainState", signature=c("ChainState","FunctionBasedBlockUpdater","BlockUpdaters","StepInfo")
    , function(chainState, blockUpdater, blockUpdaters, stepInfo) {
        #assert_that( length(getBlockIndicesIParametersToUpdate(blockUpdater))>0 )
        subSpace = getSubSpace(blockUpdater)
        iParametersToUpdate <- getBlockIndicesIParametersToUpdate(blockUpdater)
        parBounds <- getIndexedParBounds(subSpace)[iParametersToUpdate, ,drop=FALSE]
        resFUpdateBlock <- do.call( blockUpdater@fUpdateBlock, c(list(
                    getBlockParameters(chainState, blockUpdater)
                    , getBlockIndicesIParametersToUpdateInBlock(blockUpdater)
                    , lowerParBounds=parBounds[,1L]
                    , upperParBounds=parBounds[,2L]
                    , intermediates=getChainStatesIntermediates(chainState, getBlockIndicesIntermediateIdsUsed(blockUpdater))
        ), blockUpdater@argsFUpdateBlock) )
        if (resFUpdateBlock$isUpdated) {
            if(!isInIndexedParBounds(subSpace, resFUpdateBlock$xC, iParametersToUpdate)) stopDemac(
                        "Update function fUpdteBlock must return parameters within subspace."
                        ,"But was ",resFUpdateBlock$xC,", which is outside subSpace parameterBounds: lower=",getLowerParBounds(subSpace)
                        ," upper=",getUpperParBounds(subSpace)
            )
            blockParametersToUpdate(chainState,blockUpdater, blockUpdaters) <- resFUpdateBlock$xC
            isChangedByLastUpdate(chainState, blockUpdater) <- TRUE
        }
        chainState
    })



#setMethod("updateParameterBlock", signature=c("FunctionBasedBlockUpdater", "numeric")
#    , function(object, parameters, ...) {
#        assert_that( length(iParametersToUpdate(object))>0 )
#        subSpace = subSpace(object)
#        chainState = chainState(object)
#        resFUpdateBlock <- do.call( object@fUpdateBlock, c(list(parameters[iParametersUsed(object)]
#             , iParametersToUpdate(object)
#             , lowerParBounds(subSpace), upperParBounds(subSpace), chainState@intermediate)
#             , object@argsFUpdateBlock) )
#        chainState@parameters[ iParametersToUpdate(object) ] <- resFUpdateBlock$xC
#    })

#removeMethod("initialize", "FunctionBasedBlockUpdater")
#setMethod("initialize", "FunctionBasedBlockUpdater", function(.Object, ...){
#        .Object <- callNextMethod()
#        recover()
#        .Object
#    })

#library(twDev)    # automatic generation of GSetter
#---generateAndPrintS4GSetters("FunctionBasedBlockUpdater")
if(!exists("getFUpdateBlock")) setGeneric("getFUpdateBlock", function(object) standardGeneric("getFUpdateBlock"))
#' @export
setMethod("getFUpdateBlock", signature="FunctionBasedBlockUpdater", function(object) {object@fUpdateBlock})

if(!exists("getArgsFUpdateBlock")) setGeneric("getArgsFUpdateBlock", function(object) standardGeneric("getArgsFUpdateBlock"))
#' @export
setMethod("getArgsFUpdateBlock", signature="FunctionBasedBlockUpdater", function(object) {object@argsFUpdateBlock})


# do not override show: showing the function or the arguments is too long


