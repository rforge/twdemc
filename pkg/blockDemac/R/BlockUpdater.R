#' @include BlockIndices.R
#' @include SubSpace.R
#' @include ChainState.R
#' @include StepInfo.R
# need to include following statement in derived classes
# ' @include BlockUpdaters.R

#---------------------- Blockupdater class ----------------------- 
#' @export
setClass("BlockUpdater",
    contains = c("BlockIndices"), 
    representation(
        subSpace="SubSpace"                   ##<< constraints of the space of updated parameter vector
    )
    ,prototype=list(subSpace=new("SubSpace"))
)

# need to include following statement in derived classes
# ' @include BlockUpdaters.R
if(!exists("updateBlockInChainState")) setGeneric("updateBlockInChainState", function( chainState, blockUpdater, blockUpdaters, stepInfo , ...) standardGeneric("updateBlockInChainState"))

if(!exists("computeLogDensityComponents")) setGeneric("computeLogDensityComponents", function(object, ...) standardGeneric("computeLogDensityComponents"))
#' @export
setMethod("computeLogDensityComponents", signature=c("BlockUpdater")
        , function (object, parameters, intermediates, logDensityComponents, ...) {
            # needs to be implemented by subClasses
            # here as NULL object: set to -10
            logDenNames <- getLogDensityComponentNames(object)
            structure( rep(-10, length(logDenNames)), names=logDenNames)
        })



if(!exists("getLogDensityComponentNames")) setGeneric("getLogDensityComponentNames", function(object) standardGeneric("getLogDensityComponentNames"))
#' @export
setMethod("getLogDensityComponentNames", signature="BlockUpdater",
    # here return NULL-object: empty character vector. Subsclasses need to implement specific behaviour
    function(object) {
        character(0)
    }
)

if(!exists("isJumpProposalRequired")) setGeneric("isJumpProposalRequired", function(object) standardGeneric("isJumpProposalRequired"))
setMethod("isJumpProposalRequired", signature="BlockUpdater",
    # here return NULL-object: FALSE, overwrite in Metropolis-based updaters
    function(object) {
        FALSE
    }
)

if(!exists("computeChainStatesLogDensityComponents")) setGeneric("computeChainStatesLogDensityComponents", function(chainState, object, ...) standardGeneric("computeChainStatesLogDensityComponents"))
#' @export
setMethod("computeChainStatesLogDensityComponents", signature=c("ChainState","BlockUpdater")
        , function (chainState, object, ...) {
            resLogDen <- computeLogDensityComponents(object
                    , parameters= getBlockParameters(chainState, object)  
                    , intermediates=getChainStatesIntermediates(chainState, getBlockIndicesIntermediateIdsUsed(object))  
                    , logDensityComponents = getBlockLogDensityComponents(chainState, object)
                    , ...
            )
            blockLogDensityComponents(chainState,object) <- resLogDen
            chainState
        })


#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("BlockUpdater")
if(!exists("getBlockIndices")) setGeneric("getBlockIndices", function(object) standardGeneric("getBlockIndices"))
#' @export
setMethod("getBlockIndices", signature="BlockUpdater", function(object) {
            blockIndices <- copyBlockIndices( new("BlockIndices"), object )
            blockIndices
        })
if(!exists("blockIndices<-")) setGeneric("blockIndices<-", function(object,value) standardGeneric("blockIndices<-"))
setReplaceMethod("blockIndices", signature=c("BlockUpdater", "BlockIndices"), function(object, value) {
            object <- copyBlockIndices(object, value)
            object
        })

#setGeneric("getMcSetup", function(object) standardGeneric("getMcSetup"))
#setMethod("getMcSetup", signature="BlockUpdater", function(object) {object@mcSetup})
#setGeneric("mcSetup<-", function(object,value) standardGeneric("mcSetup<-"))
#setReplaceMethod("mcSetup", signature=c("BlockUpdater", "list"), function(object, value) {object@mcSetup <- value; object})

if(!exists("getSubSpace")) setGeneric("getSubSpace", function(object) standardGeneric("getSubSpace"))
#' @export
setMethod("getSubSpace", signature="BlockUpdater", function(object) {object@subSpace})
if(!exists("subSpace<-")) setGeneric("subSpace<-", function(object,value) standardGeneric("subSpace<-"))
#' @export
setReplaceMethod("subSpace", signature=c("BlockUpdater", "SubSpace"), function(object, value) {
            object@subSpace <- value; object})




