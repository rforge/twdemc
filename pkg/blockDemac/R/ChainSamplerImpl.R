#' @include ChainSampler.R
#' @include ChainState.R
#' @include StepInfoRangeImpl.R
#' @include SampleLogChainImpl.R
#' @include BlockUpdaters.R

# need to export to work on cluster
#' @export
setClass("ChainSamplerImpl"
    ,contains="ChainSampler"
    ,representation(  
        blockUpdaters = "BlockUpdaters"
        ,thin = "integer"
        ,sampleLog="SampleLogChain"
        # range information
        ,chainState="ChainState"
        ,stepInfoRange = "StepInfoRangeImpl"
        ,nInterval = "integer"
        # internally shared between private methods
        ,iGen0 = "integer" 
        ,countAcceptedInCurrentInterval = "integer"
    ),
    prototype=list(sampleLog=new("SampleLogChainImpl")),
    validity=function(object){
        if( getLength(object@blockUpdaters)==0 ) return("Need to specify blockUpdaters.")
        if( length(object@thin)!=1 ) return("Need to specify scalar thinning interval, thin.")
        return(TRUE)
    }
)

newChainSampler <- function(blockSpecifications, parameterNames, intermediateSpecifications=list()){
    if( missing(parameterNames) ) 
        parameterNames <- do.call(c,lapply(blockSpecifications, getParametersToUpdate))
    blockUpdaters <- newBlockUpdaters(blockSpecifications, intermediateSpecifications, parameterNames)
    chainSampler <- initializeBlockUpdaters(new("ChainSamplerImpl"), blockUpdaters)
    chainSampler
}

if(!exists("initializeBlockUpdaters")) setGeneric("initializeBlockUpdaters", function(object, ...) standardGeneric("initializeBlockUpdaters"))
setMethod("initializeBlockUpdaters", signature=c("ChainSamplerImpl"),
    function( object
        ,blockUpdaters
        ,thin=4L
    ) {
        object@blockUpdaters <- blockUpdaters
        object@thin <- thin
        validObject(object)
        object  
    })


setMethod("getBlockDimensions", signature="ChainSamplerImpl", function(object) { getBlockDimensions(object@blockUpdaters)})

setMethod("getIParametersThatRequireJumpProposal", signature=c("ChainSamplerImpl"),  function(object) { getIParametersThatRequireJumpProposal(object@blockUpdaters) })

if(!exists("setRangeSpecs")) setGeneric("setRangeSpecs", function(object, ...) standardGeneric("setRangeSpecs"))
setMethod("setRangeSpecs", signature=c("ChainSamplerImpl"),
    function( object
        ,chainState
        ,stepInfoRange 
        ,nInterval 
    ){
        object@chainState <- chainState  
        object@stepInfoRange <- stepInfoRange 
        object@nInterval <- nInterval
        object <- .checkRangeSpecConsistency(object)
        object <- .initializeSampleLog(object)            
        object@iGen0 <- 0L
        object            
    }
)

if(!exists(".checkRangeSpecConsistency")) setGeneric(".checkRangeSpecConsistency", function(object) standardGeneric(".checkRangeSpecConsistency"))
setMethod(".checkRangeSpecConsistency", signature=c("ChainSamplerImpl"),
    function (object) {
        blockDimensions <- getBlockDimensions(object@blockUpdaters)
        parameterNames <- getParameterNames(blockDimensions)
        logDensityComponentNames <- getLogDensityComponentNames(blockDimensions)
        diffs <-setdiff(parameterNames, names(getChainStateParameters(object@chainState)))
        if( length(diffs) ) stop("Missing parameters ind chainState:",paste(diffs,collapse=",") )
        diffs <-setdiff(logDensityComponentNames, names(getLogDensityComponents(object@chainState)) )
        if( length(diffs) ) stop("Missing logDensity components in chainState:",paste(diffs,collapse=",") )
        diffs <-setdiff(parameterNames, names(getStep(object@stepInfoRange)))
        if( object@nInterval*object@thin < getNGenerations(object@stepInfoRange))
            stop("stepInfoRange holds fewer generations than required.")
        if( object@nInterval*object@thin != getNGenerations(object@stepInfoRange))
            warning("stepInfoRange holds more generations than required.")
        if( length(logDensityComponentNames) != length(getTemperature(object@stepInfoRange)))
            stop("stepInfoRange@temperature requires nLogDensityComponents=",
                length(logDensityComponentNames),
                " entries, but has only ",
                length(getTemperature(object@stepInfoRange)))
        object
    }
)

if(!exists(".initializeSampleLog")) setGeneric(".initializeSampleLog", function(object) standardGeneric(".initializeSampleLog"))
setMethod(".initializeSampleLog", signature=c("ChainSamplerImpl"),
    function (object) {
        if( getPotentialNSample(object@sampleLog) < object@nInterval ){
            blockDim <- getBlockDimensions(object@blockUpdaters)
            object@sampleLog <- initializeLog(object@sampleLog, object@nInterval, 
                getParameterNames(blockDim), getLogDensityComponentNames(blockDim), getBlockNames(blockDim))
        } else {
            nSample(object@sampleLog) <- object@nInterval
            invalidate(object@sampleLog)
        }
        object
    }
)

if(!exists("subSpace<-")) setGeneric("subSpace<-", function(object,value) standardGeneric("subSpace<-"))
setReplaceMethod("subSpace", signature=c("ChainSamplerImpl", "SubSpace"), 
    function(object, value) {
        subSpace(object@blockUpdaters) <- value
        object
    })

#' @export
setMethod("computeChainState", signature=c("ChainSamplerImpl")
        , function (object 
                ### compute all intermediates and logDensity components for given parameter vector
                ,parameters     ##<< numceric vector of parameter values
        ) {
            chainState <- new("ChainState", parameters, intermediateIds=getIntermediateIdsUsed(object) )
            logDenNames <- getLogDensityComponentNames(getBlockDimensions(object)) 
            logDensityComponents(chainState) <- structure(rep(NA_real_, length(logDenNames)), names=logDenNames)
            chainState <- computeInvalidIntermediates(object@blockUpdaters, chainState)
            chainState <- computeInvalidChainStateDensities(object@blockUpdaters, chainState)
            chainState
        })

setMethod("sampleRange", signature=c("ChainSamplerImpl"),
    function(object) {
        object@iGen0 <- 0L
        #object@countAcceptedInInterval[] <- 0L # now in 
        object <- .initializeSampleLog(object)
        # iGen0 is increased in .computeUpdatedChainStateThinningIntervalAndCountAccepted
        while( object@iGen0 < object@nInterval*object@thin) {
            object <- .sampleInterval(object)
        }
        object
    }
)

if(!exists(".sampleInterval")) setGeneric(".sampleInterval", function(object) standardGeneric(".sampleInterval"))
setMethod(".sampleInterval", signature=c("ChainSamplerImpl"),
    function( object ) {
        chainStateGiven <- object@chainState
        object <- .computeUpdatedChainStateThinningIntervalAndCountAccepted(object)
        tmp <- object@chainState
        object@chainState <- computeInvalidChainStateDensities(object@blockUpdaters, object@chainState)
  if( any(is.na(object@chainState@logDensityComponents)))
            stop("Encountered NA-logDensities in chainState, maybe set to -Inf to signal unsuccessful calculation")                
        iInterval <- object@iGen0 %/% object@thin
        object@sampleLog <- recordSample(object@sampleLog, iSample=iInterval
            ,parameters =  getChainStateParameters(object@chainState)
            ,logDensityComponents = getLogDensityComponents(object@chainState)
            ,proportionAcceptedInInterval = object@countAcceptedInCurrentInterval / object@thin
        )
        #object <- .recordSample(object)
        object
    }
)

if(!exists(".computeUpdatedChainStateThinningIntervalAndCountAccepted")) setGeneric(".computeUpdatedChainStateThinningIntervalAndCountAccepted", function(object) standardGeneric(".computeUpdatedChainStateThinningIntervalAndCountAccepted"))
setMethod(".computeUpdatedChainStateThinningIntervalAndCountAccepted", signature=c("ChainSamplerImpl"),
    function(object) {
        chainStateGiven <- object@chainState
        updatedChainState <- object@chainState
        countAccepted <- integer( getLength(object@blockUpdaters) ) # initialize by 0
        for( iGen in object@iGen0+(1:object@thin) ){
            # tell stepInfoRange which row to return
            stepInfoRangeGenerationIndex(object@stepInfoRange) <- iGen
            # compute the update
            updatedChainState <- computeUpdatedChainState(updatedChainState, object@blockUpdaters, object@stepInfoRange)
            # increase count by 1 for accepted blocks
            countAccepted = countAccepted + 
                as.integer( .isAcceptedBlocks(object, updatedChainState))
        }
        object@chainState <- updatedChainState
        object@countAcceptedInCurrentInterval <- countAccepted
        object@iGen0 <- object@iGen0+ object@thin
        object
    }
)

if(!exists(".isAcceptedBlocks")) setGeneric(".isAcceptedBlocks", function(object, chainState) standardGeneric(".isAcceptedBlocks"))
setMethod(".isAcceptedBlocks", signature=c("ChainSamplerImpl","ChainState"),
    function( object, chainState ){
        isAcceptedBlocks <- blockUpdaterApply(object@blockUpdaters, function(blockUpdater){
                any(isChangedByLastUpdate( chainState, blockUpdater))
            })
        isAcceptedBlocks                    
    })

#.isAcceptedBlocks <- function( chainState, blockUpdaters ){
#            isAcceptedBlocks <- blockUpdaterApply(blockUpdaters, function(blockUpdater){
#                        any(isChangedByLastUpdate( chainState, blockUpdater))
#                    })
#            isAcceptedBlocks                    
#        }
#

if(!exists("getIntermediateIdsUsed")) setGeneric("getIntermediateIdsUsed", function(object,...) standardGeneric("getIntermediateIdsUsed"))
setMethod("getIntermediateIdsUsed", signature=c("ChainSamplerImpl"), 
        function(object
            ### Get the ids of all the intermediates used by the updaters.
        ) {
            getIntermediateIdsUsed(object@blockUpdaters)
        })




#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("ChainSamplerImpl")
if(!exists("getChainState")) setGeneric("getChainState", function(object) standardGeneric("getChainState"))
setMethod("getChainState", signature="ChainSamplerImpl", function(object) {object@chainState})
if(!exists("chainState<-")) setGeneric("chainState<-", function(object,value) standardGeneric("chainState<-"))
setReplaceMethod("chainState", signature=c("ChainSamplerImpl", "ChainState"), function(object, value) {object@chainState <- value; object})

if(!exists("getSampleLog")) setGeneric("getSampleLog", function(object) standardGeneric("getSampleLog"))
setMethod("getSampleLog", signature="ChainSamplerImpl", function(object) {object@sampleLog})

if(!exists("getThin")) setGeneric("getThin", function(object) standardGeneric("getThin"))
setMethod("getThin", signature="ChainSamplerImpl", function(object) {object@thin})
if(!exists("thin<-")) setGeneric("thin<-", function(object,value) standardGeneric("thin<-"))
setReplaceMethod("thin", signature=c("ChainSamplerImpl", "integer"), function(object, value) {object@thin <- value; object})



#if(!exists("getNBlock")) setGeneric("getNBlock", function(object) standardGeneric("getNBlock"))
#setMethod("getNBlock", signature="ChainSamplerImpl", function(object) {getLength(object@blockUpdaters)})
#
#if(!exists("getBlockNames")) setGeneric("getBlockNames", function(object) standardGeneric("getBlockNames"))
#setMethod("getBlockNames", signature="ChainSamplerImpl", function(object) {getNames(object@blockUpdaters)})



