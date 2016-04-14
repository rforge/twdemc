#' @include BlockIndices.R
#' @include SampleLogChain.R

#------------------ ChainState 

# ChainState is an active record, passed between processes
# make sure not to add references to objects potentially containing a lot of data

# need to export to work on cluster
#' @export
setClass("ChainState"
        ,representation(  
                parameters="numeric"
                ,isChangedByLastUpdate="logical"
                ,logDensityComponents="numeric"
                ,intermediates = "list"
                ,invalidTemplate = "list"
        )
        # default intermediates
        ,prototype = list( intermediates=list(), invalidTemplate=list() )
)


setMethod("initialize", signature=c("ChainState"),
        function( .Object
                , parameters            ##<< numeric vector (nParameter) of parameters
                , ...                   ##<< passed to callNextMethod
                , isChangedByLastUpdate ##<< logical vector (nParameter) indicating wheter parameters where changed.
                ## If its scalar, then its used for all parameters.
                ## Defaults to FALSE.
                , intermediates         ##<< named list: intermediate results
                , intermediateIds       ##<< character vector: alternative way of specifying intermediates.
        ## Entries Will be initialized to empty list.
        ){
            if( missing(parameters) || length(parameters)==0 )
                stop("You must provide named parameters for initializing ChainState")
            if( is.null(names(parameters)) )
                stop("You must provide names with the parameters for initializing ChainState.")
            nParameters <- length(parameters)
            if( missing(isChangedByLastUpdate) || length(isChangedByLastUpdate) == 0 )
                isChangedByLastUpdate <- rep(FALSE, nParameters)
            if( length(isChangedByLastUpdate) == 1 )
                isChangedByLastUpdate <- rep( isChangedByLastUpdate, nParameters)
            if( length(isChangedByLastUpdate) != nParameters )
                stop("isChangedByLastUpdate must be of same length as the parameters.")
            if( missing(intermediates) ){        
                if( missing(intermediateIds) ){
                    intermediateIds <- character(0)
                }
                intermediates <- structure( rep( list(list()), times=length(intermediateIds)), names=intermediateIds )
            }
            callNextMethod(.Object, parameters=parameters,...,
                    isChangedByLastUpdate=isChangedByLastUpdate, intermediates=intermediates)    
        })


setMethod("show", "ChainState",
        function(object){
            headParmsStr <- catNamedVector(head(object@parameters,4))
            cat("parms(",length(object@parameters),")(", headParmsStr,")",sep="")
            headLogDenStr <- catNamedVector(head(object@logDensityComponents,4))
            cat(" logDen(",length(object@logDensityComponents),")(", headLogDenStr,")",sep="")
            cat(" interm(",length(object@intermediates),")(",sep="")
            cat(head(names(object@intermediates),3),sep=",")
            cat(")\n")
        })

if(!exists("initializeFromSampleLogChain")) setGeneric("initializeFromSampleLogChain", function(object, ...) standardGeneric("initializeFromSampleLogChain"))
setMethod("initializeFromSampleLogChain", signature=c("ChainState"), 
        function(object, sampleLogChain, iSample=getNSample(sampleLogChain)) {
            object@parameters <- getParametersForSampleIndex(sampleLogChain,iSample)
            object@isChangedByLastUpdate <- logical( length(object@parameters) )
            object@logDensityComponents <- getLogDensityComponentsForSampleIndex(sampleLogChain, iSample)
            if( is.null(object@intermediates) ) object@intermediates <- list()
            object
        }
)

#if(!exists("getBlockLogDensityComponents")) setGeneric("getBlockLogDensityComponents", function(chainState, blockIndices) standardGeneric("getBlockLogDensityComponents"))
#setMethod("getBlockLogDensityComponents", signature=c("ChainState", "BlockIndices"), 
getBlockLogDensityComponents <- function(chainState, blockIndices) {
    chainState@logDensityComponents[ getBlockIndicesILogDensityComponents(blockIndices) ]
}
#)
#if(!exists("blockLogDensityComponents<-")) setGeneric("blockLogDensityComponents<-", function(object, blockIndices, value) standardGeneric("blockLogDensityComponents<-"))
#setReplaceMethod("blockLogDensityComponents", signature=c("ChainState", "BlockIndices",  "numeric"), 
"blockLogDensityComponents<-" <-   function(object, blockIndices, value) {
    object@logDensityComponents[ getBlockIndicesILogDensityComponents(blockIndices) ] <- value
    object
}
#)

#if(!exists("invalidateBlockLogDensities")) setGeneric("invalidateBlockLogDensities", function( object, blockIndices ) standardGeneric("invalidateBlockLogDensities"))
#setMethod("invalidateBlockLogDensities", signature=c("ChainState","BlockIndices"),
invalidateBlockLogDensities <-    function (object, blockIndices) {
    object@logDensityComponents[ getBlockIndicesILogDensityComponents(blockIndices) ] <- NA_real_
    object
}
#)

#if(!exists("getBlockParameters")) setGeneric("getBlockParameters", function(chainState, blockIndices) standardGeneric("getBlockParameters"))
#setMethod("getBlockParameters", signature=c("ChainState", "BlockIndices"), 
getBlockParameters <-    function(chainState, blockIndices) {
    chainState@parameters[ getBlockIndicesIParametersUsed(blockIndices) ]
}
#)

#if(!exists("blockParametersToUpdate<-")) setGeneric("blockParametersToUpdate<-", function(object, blockIndices, blockUpdaters, value) standardGeneric("blockParametersToUpdate<-"))
#setReplaceMethod("blockParametersToUpdate", signature=c("ChainState", "BlockIndices", "BlockUpdaters","numeric"), 
"blockParametersToUpdate<-" <-       function(object, blockIndices, blockUpdaters, value) {
    # Need to pass blockUpdaters object, that knows which other blocks and intermediates are to be
    # invalidated when parameters of this block is set
    #if( nargs() != 4L ) stop("Need to provide blockUpdaters.")
    #object@parameters[ getBlockIndicesIParametersToUpdate(blockIndices) ] <- value
    object@parameters[ getBlockIndicesIParametersToUpdate(blockIndices) ] <- value
    object <- invalidateBlockLogDensities(object, blockIndices)
    object <- invalidateDependents(blockUpdaters, object, blockIndices)            
    object
}
#)

#if(!exists("isChangedByLastUpdate")) setGeneric("isChangedByLastUpdate", function(chainState, blockIndices) standardGeneric("isChangedByLastUpdate"))
#setMethod("isChangedByLastUpdate", signature=c("ChainState", "BlockIndices"), 
isChangedByLastUpdate <-  function(chainState, blockIndices) {
    chainState@isChangedByLastUpdate[ getBlockIndicesIParametersToUpdate(blockIndices) ]
}
#)
#if(!exists("isChangedByLastUpdate<-")) setGeneric("isChangedByLastUpdate<-", function(chainState, blockIndices, value) standardGeneric("isChangedByLastUpdate<-"))
#setReplaceMethod("isChangedByLastUpdate", signature=c("ChainState", "BlockIndices", "logical"), 
"isChangedByLastUpdate<-" <-  function(chainState, blockIndices, value) {
    #if( length(chainState@isChangedByLastUpdate) < max(getBlockIndicesIParametersToUpdate(blockIndices)) )
    #    stop("indices do not match chainState extent")
    chainState@isChangedByLastUpdate[ getBlockIndicesIParametersToUpdate(blockIndices) ] <- value
    chainState
}
#)

#if(!exists("getChainStatesIntermediate")) setGeneric("getChainStatesIntermediate", function(object, ...) standardGeneric("getChainStatesIntermediate"))
#setMethod("getChainStatesIntermediate", signature=c("ChainState"), 
getChainStatesIntermediate <- function(object, intermediateId) {
    result <- object@intermediates[[ intermediateId ]]
    if( is.null(result) ){
        warning("intermediateId '",intermediateId,"' not known in chainState, returning empty list.")
        result <- list()
    }
    result
}
#)

#if(!exists("getChainStatesIntermediates")) setGeneric("getChainStatesIntermediates", function(object, ...) standardGeneric("getChainStatesIntermediates"))
#setMethod("getChainStatesIntermediates", signature=c("ChainState"), 
getChainStatesIntermediates <-   function(object
        ### Get the list of intermediate results for given Ids
        , intermediateIds   ##<< character vector: ids of intermediates
##<< if argument is missing, all intermediates are returned
) {
    if( missing(intermediateIds) ){
        return(object@intermediates)
    }
    result <- object@intermediates[ intermediateIds ]
    ##details<< 
    ## If some of the provided intermediateIds are not
    ## ids of chainState intermediates. A warning is issued, and 
    ## a list of only the exising intermediateIds is returned. 
    iMissing <- which( is.null(result) ) 
    if( length(iMissing) ){
        warning(pasteHead("intermediateIds '",as.vector(intermediateIds)[iMissing],"' not known in chainState, returning empty list."))
        result <- result[-iMissing]
    }
    ##details<< A named list of intermediate results
    result
}
#)

#if(!exists(".chainStatesIntermediate<-")) setGeneric(".chainStatesIntermediate<-", function(object, ..., value) standardGeneric(".chainStatesIntermediate<-"))
#setReplaceMethod(".chainStatesIntermediate", signature=c("ChainState"), 
".chainStatesIntermediate<-"  <- function(object, intermediateId, value) {
    if( nargs() != 3L ) stop("Need to provide intermediateId.")
    object@intermediates[[ intermediateId ]] <- value
    object
}
#)

#if(!exists(".invalidateChainStatesIntermediates")) setGeneric(".invalidateChainStatesIntermediates", function(object, ...) standardGeneric(".invalidateChainStatesIntermediates"))
#setMethod(".invalidateChainStatesIntermediates", signature=c("ChainState"), 
.invalidateChainStatesIntermediates <-  function(object
        ### Mark the intermediate state as being not up to date
        , intermediateIds
) {
    ##details<< It assigns an empty list. 
    ## An invalid, i.e. non-up to data,  intermediate result
    ## will be tested on by \code{!length(intermediateResults)}.
    # by making invalidTemplate an instance variable, avoid list creation
    # invalidTemplate <- list()
    for( intermediateId in intermediateIds ){
        object@intermediates[[ intermediateId ]] <- object@invalidTemplate
    }
    object
}
#)

if(!exists("getChainStateParameters")) setGeneric("getChainStateParameters", function(object) standardGeneric("getChainStateParameters"))
setMethod("getChainStateParameters", signature="ChainState", function(object) {object@parameters})

if(!exists("chainStateParameters<-")) setGeneric("chainStateParameters<-", function(object,value) standardGeneric("chainStateParameters<-"))
setReplaceMethod("chainStateParameters", signature=c("ChainState", "numeric"), 
        function(object, value) {
            if( !length(names(value)) ) names(value) <- names(object@parameters)
            object@parameters <- value
            # invalidate all blocks, i.e. logDensities, and intermediates
            object@logDensityComponents[] <- NA_real_
            object@intermediates <- lapply( object@intermediates, function(entry){ list() })
            object
            # the other way of setting parameters is "blockParametersToUpdate<-"
        })




#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("ChainState")
if(!exists("getLogDensityComponents")) setGeneric("getLogDensityComponents", function(object) standardGeneric("getLogDensityComponents"))
setMethod("getLogDensityComponents", signature="ChainState", function(object) {object@logDensityComponents})
if(!exists("logDensityComponents<-")) setGeneric("logDensityComponents<-", function(object,value) standardGeneric("logDensityComponents<-"))
setReplaceMethod("logDensityComponents", signature=c("ChainState", "numeric"), function(object, value) {object@logDensityComponents <- value; object})


if(!exists("intermediates<-")) setGeneric("intermediates<-", function(object,...,value) standardGeneric("intermediates<-"))
setReplaceMethod("intermediates", signature=c("ChainState", "list"), function(object, value) {object@intermediates <- value; object})



