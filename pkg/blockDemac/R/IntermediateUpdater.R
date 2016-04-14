#' @include IntermediateSpecification.R
#' @include DependenceIndices.R
#' @include ChainState.R

setClass("IntermediateUpdater",
        contains=c("IntermediateSpecification","DependenceIndices"),
        representation(
                intermediateId="character"
                ,requisiteUpdaterIds="character"	        ##<< names of IntermediateUpdater that this updater depends on
        )
)

if(!exists("initializeIntermediateUpdater")) setGeneric("initializeIntermediateUpdater",  function(object,...){ standardGeneric("initializeIntermediateUpdater") })
setMethod("initializeIntermediateUpdater", signature=c("IntermediateUpdater"),
        function(object, intermediateId, intermediateSpecification, parameterNames) {
            if( nargs() != 4L ) stop("Need to provide parameterNames.")
            object <- copy(object, intermediateSpecification)
            #if( !length(object@parametersUsed) ) object@parametersUsed <- parameterNames
            object@intermediateId <- intermediateId
            object@iParametersUsed <- .matchParameterNames( object@parametersUsed, parameterNames, fErrorMsgList=function(parNames){list(
                            "unknown parameters ",parNames," used by Intermediate ",intermediateId
                            )})
            object@intermediateIdsUsed <- object@intermediatesUsed # from Specification                
            object@requisiteUpdaterIds <- character()
            object
        })

if(!exists("initializeRequisiteUpdaters")) setGeneric("initializeRequisiteUpdaters",  function(object,...){ standardGeneric("initializeRequisiteUpdaters") })
setMethod("initializeRequisiteUpdaters", signature=c("IntermediateUpdater"),
        function(object, intermediateUpdaters) {
            object@requisiteUpdaterIds <- getBlockIndicesIntermediateIdsUsed(object) 
            object
        })


#if(!exists("getIntermediateFromUpdater")) setGeneric("getIntermediateFromUpdater",  function(object,...){ standardGeneric("getIntermediateFromUpdater") })
#setMethod("getIntermediateFromUpdater", signature=c("IntermediateUpdater"),
        getIntermediateFromUpdater <-  function(
                ### Get the intermediate from ChainState
                object          ##<< IntermediateUpdater object
                , chainState    ##<< numeric vector of all parameters to calibrate
                , isWarnOnEmpty=TRUE
        ){
            result <- getChainStatesIntermediate(chainState, object@intermediateId)
            if( isWarnOnEmpty && (length(result) == 0L) && !length(attr(result,"problems")) ) 
                warning("intermediate ",object@intermediateId," is not up to date.")
            result
        }
#        )

#if(!exists("letUpdaterUpdateIntermediateInChainState")) setGeneric("letUpdaterUpdateIntermediateInChainState",  function(object,...){ standardGeneric("letUpdaterUpdateIntermediateInChainState") })
#setMethod("letUpdaterUpdateIntermediateInChainState", signature=c("IntermediateUpdater"),
        letUpdaterUpdateIntermediateInChainState <-  function(
                ### Get the intermediate from ChainState
                object         ##<< IntermediateSpecification object
                ,chainState    ##<< ChainState to update
                ,updaters      ##<< named list of intermediate updaters
        ){
            # first determine which requisite intermediates are not up to date and update them
            requisiteUpdaters <- updaters[ object@requisiteUpdaterIds ]
            requisiteIntermediates <- lapply( requisiteUpdaters, getIntermediateFromUpdater, chainState=chainState, isWarnOnEmpty=FALSE)
            lengthIntermediates <- vapply(requisiteIntermediates, length, integer(1) )
            iToUpdate <- which(lengthIntermediates == 0)
            for( i in iToUpdate ){
                chainState <- letUpdaterUpdateIntermediateInChainState(requisiteUpdaters[[i]], chainState, updaters )
                requisiteIntermediates[[i]] <- getIntermediateFromUpdater(requisiteUpdaters[[i]], chainState)
            }
            parameters <- getChainStateParameters(chainState)[object@iParametersUsed]
            newIntermediate <- .letUpdaterComputeIntermediate(object, parameters, requisiteIntermediates)
            .chainStatesIntermediate(chainState, object@intermediateId) <- newIntermediate 
            chainState
        }
        #)

#if(!exists(".letUpdaterComputeIntermediate")) setGeneric(".letUpdaterComputeIntermediate",  function(object,...){ standardGeneric(".letUpdaterComputeIntermediate") })
#setMethod(".letUpdaterComputeIntermediate", signature=c("IntermediateUpdater"),
        .letUpdaterComputeIntermediate <- function(
                ### Compute the intermediate
                object          ##<< IntermediateSpecification object
                , parameters    ##<< numeric vector parameters used by intermediate
                , intermediates ##<< named list of requisite intermediate results 
        ){
            do.call( object@updateFunction
                    , c(list(parameters, intermediates)
                            ,object@argsUpdateFunction))            
        }
        #)



#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("IntermediateUpdater")
if(!exists("getIntermediateId")) setGeneric("getIntermediateId", function(object) standardGeneric("getIntermediateId"))
setMethod("getIntermediateId", signature="IntermediateUpdater", function(object) {object@intermediateId})

getIntermediatesIntermediateId <- function(object) {object@intermediateId}
