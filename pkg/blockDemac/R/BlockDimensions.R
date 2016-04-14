#' @export
setClass("BlockDimensions", contains="VIRTUAL")

if(!exists("getNBlock")) setGeneric("getNBlock", function(object) standardGeneric("getNBlock"))
if(!exists("getBlockNames")) setGeneric("getBlockNames", function(object) standardGeneric("getBlockNames"))

if(!exists("getNParameter")) setGeneric("getNParameter", function(object) standardGeneric("getNParameter"))
if(!exists("getNParameterWithProposal")) setGeneric("getNParameterWithProposal", function(object) standardGeneric("getNParameterWithProposal"))
if(!exists("getParameterNames")) setGeneric("getParameterNames", function(object) standardGeneric("getParameterNames"))
if(!exists("getIParametersUpdatedByBlock")) setGeneric("getIParametersUpdatedByBlock", function(object, blockIndex) standardGeneric("getIParametersUpdatedByBlock"))
#if(!exists("getIParametersUsedByBlock")) setGeneric("getIParametersUsedByBlock", function(object, blockIndex) standardGeneric("getIParametersUsedByBlock"))
if(!exists("isProposalRequiredForParameters")) setGeneric("isProposalRequiredForParameters", function(object, parameterIndices) standardGeneric("isProposalRequiredForParameters"))

if(!exists("getNLogDensityComponent")) setGeneric("getNLogDensityComponent", function(object) standardGeneric("getNLogDensityComponent"))
if(!exists("getLogDensityComponentNames")) setGeneric("getLogDensityComponentNames", function(object) standardGeneric("getLogDensityComponentNames"))
if(!exists("getILogDensityComponentsByBlock")) setGeneric("getILogDensityComponentsByBlock", function(object, blockIndex) standardGeneric("getILogDensityComponentsByBlock"))

if(!exists("isBlockDimensionsConsistentWith")) setGeneric("isBlockDimensionsConsistentWith", function(object, ...) standardGeneric("isBlockDimensionsConsistentWith"))
setMethod("isBlockDimensionsConsistentWith", signature="BlockDimensions", 
        function(object,blockDimensions) {
            if( !all(getParameterNames(object) == getParameterNames(blockDimensions)) ) return(FALSE)
            if( !all(getLogDensityComponentNames(object) == getLogDensityComponentNames(blockDimensions)) ) return(FALSE)
            if( !all(getBlockNames(object) == getBlockNames(blockDimensions)) ) return(FALSE)
            return(TRUE)
        })

