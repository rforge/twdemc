#' @include SampleDimensions.R
#' @include BlockDimensionsImpl.R

setClass("SampleDimensionsImpl",
    contains=c("SampleDimensions","BlockDimensionsImpl"),
    representation(
        nSamplePopulations = "integer"
        ,nChainInPopulation = "integer"
        ,nChainGroupsInPopulation = "integer"
    )
    ,prototype(nSamplePopulations=4L, nChainInPopulation=3L, nChainGroupsInPopulation=1L)
)

setMethod("show", "SampleDimensionsImpl",
    function(object){
        cat("(nPar=",getNParameter(object),", nDen=",getNLogDensityComponent(object)
        ,", nPop=",getNPopulation(object)
        ,", nChainInPop=",getNChainInPopulation(object)
        ,", nChainInGroup=",getNChainGroupsInPopulation(object),sep="")
cat(")\n")
    })

if(!exists("initializeSampleDimensionsImpl")) setGeneric("initializeSampleDimensionsImpl", function(object, ...) standardGeneric("initializeSampleDimensionsImpl"))
setMethod("initializeSampleDimensionsImpl", signature=c("SampleDimensionsImpl"), 
    function( object
        , blockDimensions
        , nSamplePopulations
        , nChainInPopulation
        , nChainGroupsInPopulation = 1L #max(1L,(object@nChainInPopulation %/% 3L))
    ){
        if( missing(blockDimensions)  )
            stop("Missing argument blockDimensions that is required to initialize SampleDimensionsImpl")
        if( missing(nSamplePopulations)  )
            stop("Missing argument nSamplePopulations (vector of length nPopulation) that required to initialize SampleDimensionsImpl")
        blockDimensions(object) <- blockDimensions
        object@nSamplePopulations <- nSamplePopulations
        if( !missing(nChainInPopulation)) object@nChainInPopulation <- nChainInPopulation
        object@nChainGroupsInPopulation <- nChainGroupsInPopulation 
        object
    })


setMethod("getNSamplePopulations", signature="SampleDimensionsImpl", function(object) {object@nSamplePopulations})
setMethod("getNPopulation", signature="SampleDimensionsImpl", function(object) {length(object@nSamplePopulations)})
setMethod("getNChain", signature="SampleDimensionsImpl", function(object) {getNPopulation(object)*object@nChainInPopulation})
setMethod("getNChainInPopulation", signature="SampleDimensionsImpl", function(object) {object@nChainInPopulation})
setMethod("getNChainGroupsInPopulation", signature="SampleDimensionsImpl", function(object) {object@nChainGroupsInPopulation})
setMethod("getIChainsForPopulation", signature=c("SampleDimensionsImpl","integer"), 
    function(object, iPopulation) {
        if( iPopulation > getNPopulation(object))
            stop("Has only ",getNPopulation(object)," populations, but tried to query population ",iPopulation)
        (iPopulation-1L)*object@nChainInPopulation + (1L:object@nChainInPopulation)        
    })
#setMethod("getIPopulationChains", signature=c("SampleDimensionsImpl"), 
#        function(object) {
#            rep(1:getNPopulation(object), each=object@nChainInPopulation)
#        })
setMethod("getIChainsForPopulations", signature=c("SampleDimensionsImpl","integer"), 
    function(object, iPopulations) {
        if( any(iPopulations > getNPopulation(object)))
            stop("Has only ",getNPopulation(object)," populations, but tried to query population ",max(iPopulations))
        do.call(c,lapply(iPopulations, function(iPopulation){
                    (iPopulation-1L)*object@nChainInPopulation + (1L:object@nChainInPopulation)            
                }))            
    })

setMethod("getIPopulationsForChains", signature=c("SampleDimensionsImpl","integer"), 
    function(object, iChain) {
        ((iChain-1L) %/% object@nChainInPopulation)+1L 
    })

if(!exists("modifySampleDimensionsImpl")) setGeneric("modifySampleDimensionsImpl", function(object, ...) standardGeneric("modifySampleDimensionsImpl"))
setMethod("modifySampleDimensionsImpl", signature=c("SampleDimensionsImpl"), 
    function( object, nSamplePopulations, nChainInPopulation, nChainGroupsInPopulation){
        #if( !missing(thin)  ) stop("setting thin from SampleDimensionsImpl is deprecated: set it in chainSampler.")#object@thin <- thin
        if( !missing(nSamplePopulations)  ) object@nSamplePopulations <- nSamplePopulations
        if( !missing(nChainInPopulation)  ) object@nChainInPopulation <- nChainInPopulation
        if( !missing(nChainGroupsInPopulation)  ) object@nChainGroupsInPopulation <- nChainGroupsInPopulation
        object
    })


#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("SampleDimensionsImpl")






