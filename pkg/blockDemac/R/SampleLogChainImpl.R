#' @include SampleLogChain.R

# need to export to work on cluster
#' @export
setClass("SampleLogChainImpl",
        contains="SampleLogChain",
        representation(  
                nSample = "integer"                    ##<< scalar integer: number of samples (matrix might be larger)
                ,parameters = "matrix"                  ##<< numeric matrix (nParm, nInterval): sampled parameters
                ,logDensityComponents = "matrix"        ##<< numeric matrix (nResComp, nInterval): result components of logDensity for those parameters
                ,proportionAcceptedInInterval = "matrix"     ##<< numeric matrix (nBlock, nInterval): sum of acceptance over all steps
        ),
        prototype=list(nSample=0L)
)

if(!exists("initializeLog")) setGeneric("initializeLog", function(object,...) standardGeneric("initializeLog"))
setMethod("initializeLog", signature="SampleLogChainImpl", 
        function(object,
                ### allocate and initialize storage (parameters, logDensityComponents, countAcceptedInterval) 
                nSample,                    ##<< scalar integer: number of samples, i.e. rows, to store
                parameterNames,              
                logDensityComponentNames,   
                blockNames
        ) {
            object@nSample <- nSample
            object@parameters <- matrix(NA_real_, ncol=nSample, 
                    nrow=length(parameterNames), 
                    dimnames=list(parameterNames,NULL))
            object@logDensityComponents <- matrix(NA_real_, ncol=nSample, 
                    nrow=length(logDensityComponentNames), 
                    dimnames=list(logDensityComponentNames,NULL))
                object@proportionAcceptedInInterval <- matrix(0, ncol=nSample,
                    nrow=length(blockNames), 
                    dimnames=list(blockNames,NULL))
                object
        }
)

if(!exists("recordSample")) setGeneric("recordSample", function(object,...) standardGeneric("recordSample"))
setMethod("recordSample", signature="SampleLogChainImpl", 
        function( object
        , iSample
        ,parameters
        ,logDensityComponents
        ,proportionAcceptedInInterval
    ){
        #assert_that(iSample <= object@nSample) # bottleneck
        object@parameters[,iSample] <- parameters
        object@logDensityComponents[,iSample] <- logDensityComponents
        object@proportionAcceptedInInterval[,iSample] <- proportionAcceptedInInterval
        object
    }
)

if(!exists("invalidate")) setGeneric("invalidate", function(object) standardGeneric("invalidate"))
setMethod("invalidate", signature="SampleLogChainImpl", function(object) {
            object@parameters[] <- NA
            object@logDensityComponents[] <- NA
            object@proportionAcceptedInInterval[] <- 0
            object
        })

if(!exists("nSample<-")) setGeneric("nSample<-", function(object,value) standardGeneric("nSample<-"))
setReplaceMethod("nSample", signature=c("SampleLogChainImpl", "integer"), 
        function(object, value) {
            if( getPotentialNSample(object) < value) 
                stop("nSample larger than space. Reinitialize before.")            
            object@nSample <- value; 
            object
        }
)

if(!exists(".applySampleLogComponents")) setGeneric(".applySampleLogComponents", function(object, FUN, ...) standardGeneric(".applySampleLogComponents"))
setMethod(".applySampleLogComponents", signature=c("SampleLogChainImpl","function"), 
        function(object, FUN, ...) {
            object@parameters <- FUN(object@parameters, ...) 
            object@logDensityComponents <- FUN(object@logDensityComponents, ...)
            object@proportionAcceptedInInterval <- FUN(object@proportionAcceptedInInterval, ...)
            object
        }
)
            

if(!exists("getPotentialNSample")) setGeneric("getPotentialNSample", function(object) standardGeneric("getPotentialNSample"))
setMethod("getPotentialNSample", signature="SampleLogChainImpl", function(object) { ncol(object@parameters) })


# TODO test
if(!exists("subSamples")) setGeneric("subSamples", function(object,...) standardGeneric("subSamples"))
setMethod("subSamples", signature="SampleLogChainImpl", function(object, iCols) {
        object@parameters <- object@parameters[,iCols ,drop=FALSE]
        object@logDensityComponents <- object@logDensityComponents[,iCols ,drop=FALSE]
        object@proportionAcceptedInInterval <- object@proportionAcceptedInInterval[,iCols ,drop=FALSE] 
        object@nSample <- length(iCols) 
        object
    })


# implementation of SampleLog
setMethod("getParameters", signature="SampleLogChainImpl", function(object) {object@parameters[,1:object@nSample, drop=FALSE]})
setMethod("getLogDensityComponents", signature="SampleLogChainImpl", function(object) {object@logDensityComponents[,1:object@nSample, drop=FALSE]})
setMethod("getProportionAcceptedInInterval", signature="SampleLogChainImpl", function(object) { object@proportionAcceptedInInterval[,1:object@nSample, drop=FALSE] })

setMethod("getParametersForSampleIndex", signature="SampleLogChainImpl", function(object,iSample) {object@parameters[,iSample]})
setMethod("getLogDensityComponentsForSampleIndex", signature="SampleLogChainImpl", function(object, iSample) {object@logDensityComponents[,iSample]})
setMethod("getProportionAcceptedInIntervalForSampleIndex", signature="SampleLogChainImpl", function(object, iSample) {object@proportionAcceptedInInterval[,iSample]})



#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("SampleLogChainImpl")
if(!exists("getNSample")) setGeneric("getNSample", function(object) standardGeneric("getNSample"))
setMethod("getNSample", signature="SampleLogChainImpl", function(object) {object@nSample})




