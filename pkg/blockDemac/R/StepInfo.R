#------------------ Interface to get information aside chainState for current step
setClass("StepInfo", contains="VIRTUAL")

if(!exists("getStep")) setGeneric("getStep", function(object,...) standardGeneric("getStep"))
if(!exists("getRExtra")) setGeneric("getRExtra", function(object,...) standardGeneric("getRExtra"))
if(!exists("getTemperature")) setGeneric("getTemperature", function(object,...) standardGeneric("getTemperature"))
if(!exists("getAcceptanceRate")) setGeneric("getAcceptanceRate", function(object,...) standardGeneric("getAcceptanceRate"))

if(!exists("getBlockTemperature")) setGeneric("getBlockTemperature", function(object, blockIndices) standardGeneric("getBlockTemperature"))
if(!exists(".blockTemperature<-")) setGeneric(".blockTemperature<-", function(object, blockIndices, value) standardGeneric(".blockTemperature<-"))

if(!exists("implementsStepInfo")) setGeneric("implementsStepInfo", function(object) standardGeneric("implementsStepInfo"))
setMethod("implementsStepInfo", signature="StepInfo", function(object) {
        blockIndices <- new("BlockIndices")
        if( !is.numeric(getBlockTemperature(object, blockIndices)) ) stop("getBlockTemperature(stepInfo,blockIndices) must return a numeric vector")
        .blockTemperature(object, blockIndices) <- numeric(0) 
        return(TRUE)
    })

