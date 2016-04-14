#' @include MetropolisBlockUpdater.R

# support for internal Metropolis decisions in fLogDensity
# see example function in InternalMetropolisBlockUpdater.R

# three addtional parameters supplied to fLogDensity function: logDensityComponentsAccepted, temperature, and isInternalRejectionUsed

# the resulting density components are checked for NA value, indicating an internal non-Accept


#' @export
setClass("InternalMetropolisBlockUpdater"
    , contains ="MetropolisBlockUpdater"	
    , representation(  
             iInternalLogDensityComponents = "integer"	##<< indices among logDensity components that are updated internally
    )
)

setMethod("computeLogDensityComponents", signature=c("InternalMetropolisBlockUpdater")
    , function (object
            ,parameters
            ,intermediates=list()
            ,logDensityComponents = rep(NA_real_, length(object@logDensityComponentNames), names=object@logDensityComponentNames)
        # reasonable defaults for no internal decision
        , temperature=rep(1, nLogDensityComponents)
        , logDensityComponentsAccepted = rep(-Inf, nLogDensityComponents)
        , isInternalRejectionUsed = FALSE
    ) {
        nLogDensityComponents <- length(logDensityComponents)
        # add three additional parameters to the function call 
        argsFLogDensity <- c(list(
                logDensityComponentsAccepted = logDensityComponentsAccepted
                ,temperature = temperature
                ,isInternalRejectionUsed = isInternalRejectionUsed
            ), getArgsFLogDensity( object))
        resLogDen <- .calcMaxConstrainedLogDensity( 
            fLogDensity = getFLogDensity(object)
            , parameters= parameters
            , intermediates = intermediates
            , logDensityComponents = logDensityComponents
            , argsFLogDensity = argsFLogDensity
            , maxLogDensity = getMaxLogDensity( object)
        )
        resLogDen
    })

setMethod("updateBlockInChainState", signature=c("ChainState","InternalMetropolisBlockUpdater","BlockUpdaters","StepInfo"),
    function(chainState, blockUpdater, blockUpdaters, stepInfo) {
        currentLogDenComp <- getBlockLogDensityComponents(chainState, blockUpdater)
        if (any(is.na(currentLogDenComp))) {
            # default useInternalRejection <- FALSE
            chainState <- computeChainStatesLogDensityComponents( chainState, blockUpdater )
        }
        proposedChainState <- .metropolisUpdatersGetProposedParametersInChainState(blockUpdater, chainState,  blockUpdaters, getStep(stepInfo))
        proposedChainState <- .computeInvalidIntermediatesForBlock( proposedChainState, blockUpdaters, blockUpdater )
        # changed invocation of  computeChainStatesLogDensityComponents
        proposedChainState <- computeChainStatesLogDensityComponents(
            proposedChainState, blockUpdater
            ,temperature = getBlockTemperature(stepInfo, blockUpdater)
            ,logDensityComponentsAccepted = currentLogDenComp
            ,isInternalRejectionUsed = TRUE
        )
        # NA in proposedLogDenComp in result of computeChainStatesLogDensityComponents indicates rejection
        proposedLogDenComp <- getBlockLogDensityComponents(proposedChainState, blockUpdater)
        isAccepted <- !any(is.na(proposedLogDenComp)) && 
             .metropolisDecision2(
                logDenCompCurrent = currentLogDenComp
                ,logDenCompProposal = proposedLogDenComp
                ,rExtra = getRExtra(stepInfo)
                ,tempResCompC = getBlockTemperature(stepInfo, blockUpdater)
                ,iInternalLogDensityComponents = getIInternalLogDensityComponents(blockUpdater)
            )
        currentChainState <- if( isAccepted ){
                isChangedByLastUpdate(proposedChainState, blockUpdater) <- TRUE
                proposedChainState
            } else {
                isChangedByLastUpdate(chainState, blockUpdater) <- FALSE
                chainState
            }
        currentChainState
    })


#library(twDev)    # automatic generation of GSetter
#---generateAndPrintS4GSetters("InternalMetropolisBlockUpdater")
if(!exists("getIInternalLogDensityComponents")) setGeneric("getIInternalLogDensityComponents", function(object) standardGeneric("getIInternalLogDensityComponents"))
setMethod("getIInternalLogDensityComponents", signature="MetropolisBlockUpdater", function(object) {object@iInternalLogDensityComponents})


    


