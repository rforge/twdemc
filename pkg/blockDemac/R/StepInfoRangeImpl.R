#' @include StepInfo.R

# need to export to work on cluster
#' @export
setClass("StepInfoRangeImpl",
    contains="StepInfo",	    # implements interface
    representation(  
        iGeneration="integer"	        # local counter: the generation for which to return rows of the stepInfo
        ,step="matrix"
        ,rExtra="numeric"
        ,temperature = "matrix"
        ,acceptanceRate = "numeric"	    # numeric scalar: current acceptance rate: same for entire range
    )
    ,prototype( acceptanceRate=0.25)
)

setMethod("getStep", signature="StepInfoRangeImpl", function(object) {object@step[,object@iGeneration]})
setMethod("getRExtra", signature="StepInfoRangeImpl", function(object) {object@rExtra[object@iGeneration]})
setMethod("getTemperature", signature="StepInfoRangeImpl", function(object) {object@temperature[,object@iGeneration]})
setMethod("getAcceptanceRate", signature="StepInfoRangeImpl", function(object) {object@acceptanceRate})


setMethod("initialize", signature=c("StepInfoRangeImpl"),
    function( .Object
        , step
        , rExtra
        , temperature
        , nLogDensityComponents
        , nGeneration
        , iGeneration
        , ...
    ){
        # allow emty constructor for subclasses
        if( missing(step) && missing(rExtra) && missing(temperature) && missing(nLogDensityComponents) && missing(nGeneration) && missing(iGeneration) ){
            return( callNextMethod(.Object,  ...) )
        }
        # if only a vector is given as step, then repeat it nGeneration
        if( !missing(step) && !is.matrix(step) ) step <- matrix(step, 
                    nrow=length(step)
                    , ncol={if(missing(nGeneration)) 1 else nGeneration}
                    , dimnames=list(names(step),NULL))
        if( missing(step) || !is.matrix(step) || nrow(step)==0 )
            stop("Need to provide a parameter step as matrix (nParm, nGen) when constructing a StepInfoRangeImpl.")
        if( missing(nGeneration)) nGeneration <- ncol(step)
        if( ncol(step) != nGeneration ) stop("number of columns in step must match nGeneration")
        if( !length(rownames(step)) ) stop("must provide (row)names for parameters")            
        #
        if( missing(rExtra) || length(rExtra)==0 )
            rExtra <- rep(0, nGeneration)
        #
        if( !missing(temperature) && !is.matrix(temperature)) temperature <- matrix(temperature
                  , nrow={if( missing(nLogDensityComponents)) length(temperature) else nLogDensityComponents}
                  , ncol=nGeneration)
        # if only a vector is given, assume temperature not to change between generations
        if( missing(temperature) || length(temperature)==0 ) temperature <- matrix(1
                    , nrow={if( missing(nLogDensityComponents)) 1 else nLogDensityComponents}
                    , ncol=nGeneration)
        if( missing(nLogDensityComponents) ) nLogDensityComponents <- nrow(temperature)
        if( nrow(temperature) != nLogDensityComponents ) stop("number of rows in temperature must correspond to nLogDensityComponents.")
        if( ncol(temperature) != nGeneration ) stop("number of columns in temperature must match nGeneration")
        #
        if( missing(iGeneration) ) iGeneration <- 1L
        callNextMethod(.Object, step=step, rExtra=rExtra, temperature = temperature, iGeneration=iGeneration, ...)    
    })

if( !exists("getNGenerations")) setGeneric("getNGenerations", function(object) standardGeneric("getNGenerations"))
setMethod("getNGenerations", signature=c("StepInfoRangeImpl"), function(object) { ncol(object@step) })
            


setMethod("show", "StepInfoRangeImpl",
    function(object){
        if( nrow(object@step) == 0 ) cat("not initialized\n") else {
            stepParmsStr <- catNamedVector(head(object@step[,1],4))
            cat("step(",paste(dim(object@step),collapse=','),")(", stepParmsStr,")",sep="")
            cat(" rExtra=",head(object@rExtra,4),sep="")
            headTemperatureStr <- paste(head(object@temperature[1,],4),collapse=",")
            cat(" temp(",paste(dim(object@temperature),collapse=','),")(", headTemperatureStr,")",sep="")
            cat("\n")
        }
    })

setMethod("getBlockTemperature", signature=c("StepInfoRangeImpl", "BlockIndices"), 
    function(object, blockIndices) {
        object@temperature[getBlockIndicesILogDensityComponents(blockIndices), object@iGeneration]
    }
)
setReplaceMethod(".blockTemperature", signature=c("StepInfoRangeImpl", "BlockIndices",  "numeric"), 
    function(object, blockIndices, value) {
        object@temperature[getBlockIndicesILogDensityComponents(blockIndices), object@iGeneration ] <- value
        object
    }
)

if(!exists("rangeTemperature<-")) setGeneric("rangeTemperature<-", function(object,value) standardGeneric("rangeTemperature<-"))
setReplaceMethod("rangeTemperature", signature=c("StepInfoRangeImpl", "matrix"), function(object, value) {object@temperature <- value; object})

if(!exists("rangeAcceptanceRate<-")) setGeneric("rangeAcceptanceRate<-", function(object,value) standardGeneric("rangeAcceptanceRate<-"))
setReplaceMethod("rangeAcceptanceRate", signature=c("StepInfoRangeImpl", "numeric"), function(object, value) {object@acceptanceRate <- value; object})


#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("StepInfoRangeImpl")

if(!exists("getGenerationIndex")) setGeneric("getGenerationIndex", function(object) standardGeneric("getGenerationIndex"))
setMethod("getGenerationIndex", signature="StepInfoRangeImpl", function(object) {object@iGeneration})

if(!exists("generationIndex<-")) setGeneric("generationIndex<-", function(object,value) standardGeneric("generationIndex<-"))
setReplaceMethod("generationIndex", signature=c("StepInfoRangeImpl", "integer"), function(object, value) {object@iGeneration <- value; object})

# add non-dispatching function calls, because called for each generation
"stepInfoRangeGenerationIndex<-" <- function(object, value) {object@iGeneration <- value; object}



    
