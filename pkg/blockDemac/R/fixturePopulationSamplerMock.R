#' @include PopulationSampler.R
#' @include sampleRange.R

# initialize sample after setting initial state (i.e. chainStates) 
# and return portions of it during sampleRangeAllChains

setClass("PopulationSamplerMock"
    , contains="PopulationSampler"
    ,representation( 
        sampleLogChains = "list"
    )
)

setMethod(".samplePopulations", signature=c("PopulationSamplerMock"),
    function( object, ...
    ){
        # now chainStates are initialized
        object <- initializeMockSample(object)
        object <- callNextMethod(object, ...) 
        object
    })

if(!exists("initializeMockSample")) setGeneric("initializeMockSample", function(object,...) standardGeneric("initializeMockSample"))
setMethod("initializeMockSample", signature="PopulationSamplerMock", 
    function(object,...){
        sampleDimensions <- getSampleDimensions(object@sampleLogs)
        x0Chains <- abind( lapply( object@chainStates, getChainStateParameters ), along=2 )
        logDensityComponents0 <- abind( lapply( object@chainStates, getLogDensityComponents ), along=2)
        if( length(logDensityComponents0) && all(is.na(logDensityComponents0)) )
            logDensityComponents0 <- structure(-10-(1:getNLogDensityComponent(sampleDimensions)), names=getLogDensityComponentNames(sampleDimensions)) 
        fx <- .fixtureSampleLogs(
            sampleDimensions = sampleDimensions
            ,x0 = x0Chains
            ,logDensityComponents0 = logDensityComponents0
        )        
        object@sampleLogChains <- do.call(c, sapply(1:getNPopulation(sampleDimensions), function(iPop){
                    getSampleLogChainsForPopulations(fx$sampleLogs,iPop)
                }))
        object
    })

#------------------------ 
# sampleLog will be x0*iStep
# countAcceptedInInterval=1L for all 
setMethod("sampleRangeAllChains", signature="PopulationSamplerMock", 
    function( object                
        ,chainStates                         ##<< list with current state of each chain taking part in step (changes between generations)
        ,intervalInfoChains = list(         ##<< list with information on interval     (changes between nIntervals)
            ##describe<<
            ,iChainsStep = integer(0)   ##<< integer vector (nChainRange) of indices of chains that take part in this step, that each chain belongs to.
            ,iSample0=0L                ##<< +1L gives the current sample, i.e. thinning interval (used to index temperature)
            ,nSample=3L                 ##<< number of samples to generate per chain in this range
            ,step = numeric(0)           ##<< numeric array (nPar, nInterval*thin, nChainRange) of poposed jumps (as differences)
            ,rExtra = numeric(0)            ##<< numeric matrix (nInterval*thin, nChainRange): Snooker's update Metropolis rule correction  
            ,pAccept = numeric(0)           ##<< vector (nBlock, nChainRange): current acceptance rate
        ##end<<
        )
        ,mcCommon=list() ##<< with entry sampleDimensions
    ) {
        thin <- getThin(mcCommon$chainSampler)
        #iSample0 <- ((intervalInfoChains$iGenerations[1]-1) %/% thin) 
        #nSample <- length(intervalInfoChains$iGenerations) %/% thin
        iSample0 <- intervalInfoChains$iSample0
        nSample <- intervalInfoChains$nSample
        sampleLogChains <- object@sampleLogChains[intervalInfoChains$iChainsStep]
        sampleLogChainsInterval <- lapply( sampleLogChains, subSamples, iSample0+(1:nSample))
        lapply( seq_along(intervalInfoChains$iChainsStep), function(i){
                list(
                    chainState = initializeFromSampleLogChain( chainStates[[i]], sampleLogChainsInterval[[i]] )
                    ,sampleLog = sampleLogChainsInterval[[i]]
                )
            })
        #value<<  list for each chain, with each entry a result of \code{\link{.sampleRangeOneChain}}
    }
)
