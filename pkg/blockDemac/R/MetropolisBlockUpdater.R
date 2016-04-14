#' @include BlockUpdater.R
#' @include BlockUpdaters.R

#' @export
setClass("MetropolisBlockUpdater"
        , contains ="BlockUpdater"	    # chainStateIndices, subSpace
        , representation(  
                fLogDensity="function"	    ##<< f(theta, iParms, upperParBounds, lowerParBounds, intermediate, ...) 
                , argsFLogDensity = "list"	##<< further arguments to fLogDensity
                , logDensityComponentNames = "character"	##<< names of the vector result of fLogDen
                , maxLogDensity = "numeric"	##<< numeric vector (nResComp): maximum logDensity (overfitting control, usually -1/2 nObs). Order is important
                , isMarkUpdatedWhenNotAccepted = "logical"	##<< set to TRUE to mark as updated, even if not accepted (to trigger recomputation involving some randomness)
        )
        ,validity=function(object){
            if( length(object@logDensityComponentNames)==0 ) return("Need to provide Names of LogDensity components.")
            if( is.null(body(object@fLogDensity)) ) return("Need to provide argument fLogDensity, a function that calculates the log-Density.")
            if( length(object@maxLogDensity) != length(object@logDensityComponentNames) ) return("length of maxLogDensity must match number of logDensity components.")
            return(TRUE)
        }
        ,prototype(isMarkUpdatedWhenNotAccepted=FALSE)
)

#' @export
setMethod("getLogDensityComponentNames", signature="MetropolisBlockUpdater", function(object) {object@logDensityComponentNames})

# requires jumps to be proposed for parameters updated by this block
setMethod("isJumpProposalRequired", signature="MetropolisBlockUpdater",  function(object) {  TRUE  })

#' @export
setMethod("computeLogDensityComponents", signature=c("MetropolisBlockUpdater")
        , function (object 
                ,parameters
                ,intermediates=list()
                ,logDensityComponents = rep(NA_real_, length(object@logDensityComponentNames), names=object@logDensityComponentNames)
        ) {
            .metropolisUpdatersComputeLogDensityComponents(object, parameters, intermediates, logDensityComponents )
        })

.metropolisUpdatersComputeLogDensityComponents <- function (object 
                ,parameters
                ,intermediates=list()
                ,logDensityComponents = rep(NA_real_, length(object@logDensityComponentNames), names=object@logDensityComponentNames)
        ) {
            resLogDen <- .calcMaxConstrainedLogDensity( 
                    fLogDensity = object@fLogDensity
                    , parameters = parameters
                    , intermediates = intermediates
                    , logDensityComponents = logDensityComponents  
                    , argsFLogDensity = object@argsFLogDensity
                    , maxLogDensity = object@maxLogDensity
            )
            resLogDen
        }

#removeMethod("initialize", signature=c("MetropolisBlockUpdater"))
setMethod("initialize", signature=c("MetropolisBlockUpdater"),
        function( .Object, ...,maxLogDensity, logDensityComponentNames){
            # allow calling without arguments to allow subclasses
            if( missing(logDensityComponentNames) && missing(maxLogDensity) && length(list(...))==0)
                return(callNextMethod(.Object))
            if( missing(maxLogDensity) || length(maxLogDensity)==0 ){
                if( missing(logDensityComponentNames) ) stop("must specify logDensityComponentNames when creating MetropolisBlockUpdater")
                # do not do overfitting control, without explicitely requiring it                    
                maxLogDensity <- rep(+Inf, length(logDensityComponentNames)) 
                #maxLogDensity <- numeric( length(logDensityComponentNames) )
            }
            .Object <- callNextMethod(.Object, ...,maxLogDensity=maxLogDensity, logDensityComponentNames=logDensityComponentNames)
        })


#' @export
setMethod("updateBlockInChainState", signature=c("ChainState","MetropolisBlockUpdater","BlockUpdaters","StepInfo"),
        function(chainState, blockUpdater, blockUpdaters, stepInfo) {
            #assert_that( length(getBlockIndicesIParametersToUpdate(blockUpdater))>0 )
            chainStateGiven <- chainState
            currentLogDenComp <- getBlockLogDensityComponents(chainState, blockUpdater)
            if (any(is.na(currentLogDenComp))) {
                # need to stay S4, because required from all BlockUpdaters
                chainState <- .metropolisUpdatersComputeChainStatesLogDensityComponents(chainState, blockUpdater)
            }
            proposedChainState <- .metropolisUpdatersGetProposedParametersInChainState(blockUpdater, chainState, blockUpdaters, getStep(stepInfo))
            isAccepted <- tryCatch({
                        # assigning new parameters will invalidate intermediates and blocks that depend on
                        proposedChainState <- .computeInvalidIntermediatesForBlock( proposedChainState, blockUpdaters, blockUpdater )
                        # when computing intermediates, may use stopDemacInvalidChainState() 
                        # to signal non-acceptance instead of stopping program, 
                        # e.g. in ioSpec of GP compileDiscrepancyBlocksForStream:
                        # stopDemacInvalidChainState("iOSpec: given logPsi yield infinite psi.")
                        proposedChainState <- .metropolisUpdatersComputeChainStatesLogDensityComponents(proposedChainState, blockUpdater)
                        #if( names(blockUpdater@iParametersToUpdate)[1] == "logPsi_rich") recover()            
                        logDenCompCurrent = getBlockLogDensityComponents(chainState, blockUpdater)
                        logDenCompProposal = getBlockLogDensityComponents(proposedChainState, blockUpdater)
                        isAccepted <- .metropolisDecision2(
                                logDenCompCurrent = logDenCompCurrent
                                ,logDenCompProposal = logDenCompProposal
                                ,rExtra = getRExtra(stepInfo)
                                ,tempResCompC = getBlockTemperature(stepInfo, blockUpdater)
                                ,iInternalLogDensityComponents = integer(0)
                        )
                        isAccepted
                    }
                    , demacInvalidChainState=function(e){
                        warning(e$message)
                        FALSE   # isAccepted is FALSE, but do not break
                    } 
                   )
            currentChainState <- if( isAccepted ){
    #print(logDenCompCurrent); print(logDenCompProposal)
    #recover()         
                logDenCompProp <- getBlockLogDensityComponents(proposedChainState,blockUpdater) 
                if( (any(is.na(logDenCompProp))) ) stop("MetropolisBlockUpdater: encountered NA logDenCompProp.") 
                isChangedByLastUpdate(proposedChainState, blockUpdater) <- TRUE
                proposedChainState
            } else {
                isChangedByLastUpdate(chainState, blockUpdater) <- blockUpdater@isMarkUpdatedWhenNotAccepted
                chainState
            }
            currentChainState
        })

# copied from S4 method in BlockUpdater, to avoid dispatching overhead
.metropolisUpdatersComputeChainStatesLogDensityComponents  <- function (chainState, object, ...) {
            resLogDen <- .metropolisUpdatersComputeLogDensityComponents(object
                    , parameters= getBlockParameters(chainState, object)  
                    , intermediates=getChainStatesIntermediates(chainState, getBlockIndicesIntermediateIdsUsed(object))  
                    , logDensityComponents = getBlockLogDensityComponents(chainState, object)
                    , ...
            )
            blockLogDensityComponents(chainState,object) <- resLogDen
            chainState
        }

        
#if(!exists("getProposedParametersInChainState")) setGeneric("getProposedParametersInChainState", function(object,chainState,blockUpdaters,step) standardGeneric("getProposedParametersInChainState"))
#setMethod("getProposedParametersInChainState", signature=c("MetropolisBlockUpdater","ChainState","BlockUpdaters","numeric"),
        .metropolisUpdatersGetProposedParametersInChainState <-  function(object,chainState,blockUpdaters,step){
            #assert_that ( length(step)==length( getChainStateParameters(chainState)) )
            subSpace = getSubSpace(object)
            iParametersToUpdate <- getBlockIndicesIParametersToUpdate(object) 
            newBlockParameters <- newBlockParameters0 <- getChainStateParameters(chainState)[iParametersToUpdate] + step[iParametersToUpdate]
            newBlockParameters <- reflectAtIndexedBounds(subSpace, newBlockParameters0, iParametersToUpdate )
            #if( (names(newBlockParameters)[1]=="logPsi_rich") && (newBlockParameters < -2.82 )) recover()            
            # measured performance, but calling non S4 dispatch of reflectBoundaries is not better
            #newBlockParameters <- .letSubSpaceReflectBoundaries(subSpace, newBlockParameters )
            # MAYBE apply discretization function
            # assigning block parameters will invalidate depending intermediates and blocks
            blockParametersToUpdate(chainState, object, blockUpdaters) <- newBlockParameters         
            chainState
        }
        



#library(twDev)    # automatic generation of GSetter
#---generateAndPrintS4GSetters("MetropolisBlockUpdater")
if(!exists("getFLogDensity")) setGeneric("getFLogDensity", function(object) standardGeneric("getFLogDensity"))
#' @export
setMethod("getFLogDensity", signature="MetropolisBlockUpdater", function(object) {object@fLogDensity})

if(!exists("getArgsFLogDensity")) setGeneric("getArgsFLogDensity", function(object) standardGeneric("getArgsFLogDensity"))
#' @export
setMethod("getArgsFLogDensity", signature="MetropolisBlockUpdater", function(object) {object@argsFLogDensity})

if(!exists("getMaxLogDensity")) setGeneric("getMaxLogDensity", function(object) standardGeneric("getMaxLogDensity"))
#' @export
setMethod("getMaxLogDensity", signature="MetropolisBlockUpdater", function(object) {object@maxLogDensity})




