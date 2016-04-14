if(!exists("sampleRangeAllChains")) setGeneric("sampleRangeAllChains", function(object,...) standardGeneric("sampleRangeAllChains"))
setMethod("sampleRangeAllChains", signature="PopulationSampler", 
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
                ,mcCommon
                ,cluster=object@cluster                 ##<< a cluster object, e.g. \code{\link{makeCluster}(detectCores())}
        ) {
            nChain <- length(intervalInfoChains$iChainsStep)
            mcCommonList <- if (.isParallel(cluster)) list() else list(mcCommon=mcCommon)
            iPopulationsStep <- getIPopulationsForChains(mcCommon$sampleDimension, intervalInfoChains$iChainsStep)
            argsChains <- lapply(1:nChain, function(iChain){
                        iPopulation <- iPopulationsStep[iChain] 
                        c( mcCommonList, list(
                                        chainState = chainStates[[iChain]]
                                        ,iPopulation = iPopulation
                                        ,iSample0 = intervalInfoChains$iSample0
                                        ,nSample=intervalInfoChains$nSample
                                        ,stepInfoRange = new("StepInfoRangeImpl"
                                                ,step = adrop(intervalInfoChains$step[,,iChain ,drop=FALSE],drop=3L)                            
                                                ,rExtra = intervalInfoChains$rExtra[,iChain]
                                                ,acceptanceRate = intervalInfoChains$pAccept[iChain]
                                                #,temperature = mcCommon$temperatures[[ iPopulation ]][,
                                                #        intervalInfoChains$iGenerations ,drop=FALSE]                           
                                                # temperature is set in .sampleRangeOneChain (to avoid passing too much data)
                                        )
                                ))
                    })
            #argsChain <- argsChains[[1]]
            doCallSampleRangeOneChain <- function(argsChain){
                tmp <- do.call(.sampleRangeOneChain, argsChain)
            }
            # for debugging errors in cost function in parallel, initialize updater with: fLogDensity=dumpOnError(<origFLogDensity>) 
            # for debugging blockDemac remote
            #res <- parLapplySave(cluster, argsChains, dumpOnError(doCallSampleRangeOneChain))
            res <- parLapplySave(cluster, argsChains, doCallSampleRangeOneChain)
            #value<<  list for each chain, with each entry a result of \code{\link{.sampleRangeOneChain}}
        }
)

# need to be public (exported from package) in order to work on cluster
#' @export
.sampleRangeOneChain <- function(
        ### sampling of a range of samples from one chain, done on one node
        chainState,     ##<< ChainState: current state of the chain 
        iPopulation,    ##<< scalar integer: index of the population that the chain belongs to
        iSample0,
        nSample,
        stepInfoRange,  ##<< StepInfoRange: proposed jumps and rExtra for all generations
        mcCommon = list() ##<< list: information that does not change during MC run: \itemize{
## \item SampleDimensions: SampleDimensions object, configured with blocks and thin
## \item chainSampler: ChainSamplerImpl object
## \item subSpaces: list(nPop) with SubSpace objects
## \item temperatures: list(nPop) with temperature for all generations:  matrix (nGeneration, nLogDensityComponent)
## }
) {
    ##details<< 
    ## Argument \code{mcCommon} may hold large amounts of data, 
    ## e.g. as arguments to the block updaters.
    ## It does not need to be transferred to cluster nodes each call of this
    ## functions. Rather, if it is exported once by 
    ## \code{sfExport("mcCommon", ".sampleRangeOneChain", local=FALSE)} or
    ## \code{clusterExport(cl, c("mcCommon", ".sampleRangeOneChain"), environment())} 
    ## and retrieved from globalenv in this function. 
    if( !length(mcCommon) ){
        recover()
        mcCommon <- get("mcCommon", envir=globalenv() )  
    }
    assert_that( all(c("chainSampler","subSpacePopulations","temperatures","sampleDimensions") %in% names(mcCommon)))
    sampleDimensions <- mcCommon$sampleDimensions    # with blocks and thin set
    chainSampler <- mcCommon$chainSampler    # with blocks and thin set
    #nInterval <- length(iGenerations) %/% getThin(sampleDimensions)
    subSpace(chainSampler) <- mcCommon$subSpacePopulations[[iPopulation]]
    thin <- getThin(chainSampler)
    iGenerations = (iSample0)*thin+(1:(nSample*thin))
    rangeTemperature(stepInfoRange) <- mcCommon$temperatures[[iPopulation]][,iGenerations ,drop=FALSE]
    chainSampler <-  setRangeSpecs(chainSampler,
            chainState=chainState, nInterval=nSample, stepInfoRange=stepInfoRange)
    updatedChainSampler <- sampleRange(chainSampler)
    list(
            chainState = getChainState(updatedChainSampler)
            ,sampleLog = getSampleLog(updatedChainSampler)
    )
}

