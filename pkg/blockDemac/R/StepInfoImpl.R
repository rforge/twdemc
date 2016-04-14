#' @include StepInfo.R

setClass("StepInfoImpl",
    contains="StepInfo",	    # implements interface
    representation(  
        step="numeric"
        ,rExtra="numeric"
        ,temperature = "numeric"
        ,acceptanceRate = "numeric"
    )
)

setMethod("getStep", signature="StepInfoImpl", function(object) {object@step})
setMethod("getRExtra", signature="StepInfoImpl", function(object) {object@rExtra})
setMethod("getTemperature", signature="StepInfoImpl", function(object) {object@temperature})
setMethod("getAcceptanceRate", signature="StepInfoImpl", function(object) {object@acceptanceRate})


setMethod("initialize", signature=c("StepInfoImpl"),
    function( .Object, step, rExtra, temperature, nLogDensityComponents, ...){
        # allow emty constructor for subclasses
        if( missing(step) && missing(rExtra) && missing(temperature) && missing(nLogDensityComponents)  ){
            return( callNextMethod(.Object,  ...) )
        }
        if( missing(step) || length(step)==0 )
            stop("Need to provide parameter step when constructing a StepInfoImpl.")
        if( missing(rExtra) || length(rExtra)==0 )
            rExtra <- 0
        if( missing(temperature) || length(temperature)==0 )
            temperature <- 1
        if( length(temperature)==1 ){
            if( missing(nLogDensityComponents) ) nLogDensityComponents <- 1
            temperature <- rep(temperature, nLogDensityComponents )
        }
        callNextMethod(.Object, step=step, rExtra=rExtra, temperature = temperature, ...)    
    })

    


setMethod("show", "StepInfoImpl",
    function(object){
        stepParmsStr <- catNamedVector(head(object@step,4))
        cat("step(",length(object@step),")(", stepParmsStr,")",sep="")
        cat(" rExtra=",object@rExtra,sep="")
        headTemperatureStr <- catNamedVector(head(object@temperature,4))
        cat(" temp(",length(object@temperature),")(", headTemperatureStr,")",sep="")
        cat("\n")
    })

setMethod("getBlockTemperature", signature=c("StepInfoImpl", "BlockIndices"), 
    function(object, blockIndices) {
        object@temperature[ getBlockIndicesILogDensityComponents(blockIndices)  ]
    }
)
setReplaceMethod(".blockTemperature", signature=c("StepInfoImpl", "BlockIndices",  "numeric"), 
    function(object, blockIndices, value) {
        object@temperature[ getBlockIndicesILogDensityComponents(blockIndices) ] <- value
        object
    }
)

    
