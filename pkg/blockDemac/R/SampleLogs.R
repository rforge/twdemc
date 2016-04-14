#' @include SampleLogPopulation.R
#' @include SampleDimensionsImpl.R

#' @export
setClass("SampleLogs",
    contains="SampleDimensions",
    representation(
        sampleDimensions = "SampleDimensions"    
        ,initialPopulation = "list"	    ##<< list (nPopulation) of numeric matrices (nParm, nStepInitial, nChain) of samples in initial population 
        ,populationLogs = "list"	    ##<< list (nPopulation) of SampleLogPopulation
    #,nSample="integer"
    ),
    prototype=list(sampleDimensions=new("SampleDimensionsImpl"))
)

setMethod("initialize", signature=c("SampleLogs"),
    function( .Object
        , initialPopulation
        , ...
    ){
        .Object <- callNextMethod(.Object, ...)
        if( !missing(initialPopulation) ){
            initialPopulation(.Object) <- initialPopulation  
        }
        if( "sampleDimensions" %in% names(list(...)) && !("populationLogs" %in% names(list(...))) )
            .Object <- initializeSampleLogs(.Object)
        .Object 
    })

if(!exists("initialPopulation<-")) setGeneric("initialPopulation<-", function(object,value) standardGeneric("initialPopulation<-"))
setReplaceMethod("initialPopulation", signature=c("SampleLogs", "ANY"), 
    function(object, value) {
        object@initialPopulation <- .getCompletedInitialPopulation(object, value)
        if( .isPopulationLogsInitialized(object) ){
            # replace parameters in populations
            #stop("setting initial population on already initialized logs not implemented yet.")
            object@populationLogs <- mapply( function(populationLog, initialPopulationPop){
                    initialPopulation(populationLog) <- initialPopulationPop
                    populationLog
                },object@populationLogs, object@initialPopulation, SIMPLIFY=FALSE)
        }
        object
    })

if(!exists(".getCompletedInitialPopulation")) setGeneric(".getCompletedInitialPopulation", function(object,...) standardGeneric(".getCompletedInitialPopulation"))
setMethod(".getCompletedInitialPopulation", signature=c("SampleLogs"), 
    function( object, initialPopulation ){
        nPop <- getNPopulation(object@sampleDimensions)
        if( nPop == 0) stop("Need to set sampleDimensions before setting initialPopulation.")
        nChainInPop <- getNChainInPopulation(object@sampleDimensions)
        nParm <- getNParameter(object@sampleDimensions)
        parNames <- getParameterNames(object@sampleDimensions)
        if( is.array(initialPopulation) ){
            ##details<<
            ## Initial population can be specified as one array across all chains 
            ## instead of a list of matrices for each chain
            initialPopulation <- lapply( 1:nPop, function(iPop){
                    iChains <- getIChainsForPopulation(object@sampleDimensions, iPop)
                    initialPopulation[,, iChains,drop=FALSE]
                })
        }
        if( missing(initialPopulation) || !length(initialPopulation) ){ 
            initialPopulation <- lapply(1:nPop, function(iPop){
                    array( NA_real_, dim=c(nParm, 0, nChainInPop), dimnames=list(parNames,NULL,NULL))
                })
        }
        if( length(initialPopulation) != nPop )
            stop("An initial population must be provided for each population")
        for( iPop in 1:nPop ){ 
            if( nrow(initialPopulation[[iPop]]) != nParm ) stop("mismatch in number of parameters in initial population.")
            rownames(initialPopulation[[iPop]]) <- parNames
        }
        initialPopulation
    }
)

if(!exists(".isPopulationLogsInitialized")) setGeneric(".isPopulationLogsInitialized", function(object,...) standardGeneric(".isPopulationLogsInitialized"))
setMethod(".isPopulationLogsInitialized", signature=c("SampleLogs"), function(object){ length(object@populationLogs) != 0L })

setMethod("show", "SampleLogs",
    function(object){
        show(object@sampleDimensions)
        cat("\n")
    })

if(!exists("sampleDimensions<-")) setGeneric("sampleDimensions<-", function(object,value) standardGeneric("sampleDimensions<-"))
setReplaceMethod("sampleDimensions", signature=c("SampleLogs", "SampleDimensions"), 
    function(object, value) {
        #if( any(getNSamplePopulations(value)==0) ) stop("Zero record sampling log not supported.")
        object@sampleDimensions <- value
        object@initialPopulation <- list()
        object <- initializeSampleLogs(object)
        object
    })

if(!exists("initializeSampleLogs")) setGeneric("initializeSampleLogs", function(object,...) standardGeneric("initializeSampleLogs"))
setMethod("initializeSampleLogs", signature="SampleLogs", 
    function( object){
        nSamplePopulations <- getNSamplePopulations(object@sampleDimensions)
        object@populationLogs <- list() # marks as uninitialized
        if( !length(nSamplePopulations) ) stop("Need to set sampleDimensions before initializing SampleLogs.")
        # if initialPopulation has not been set yet, initialize to zero record initialPopulation 
        if( !length(object@initialPopulation) ) initialPopulation(object) <- numeric(0)
        nChainInPop <- getNChainInPopulation(object@sampleDimensions)
        object@populationLogs <- lapply( 1:getNPopulation(object), function(iPop){
                initializeLog( sampleLog <- new("SampleLogPopulation")
                    ,nSample=nSamplePopulations[iPop]
                    ,parameterNames = getParameterNames(object)
                    ,logDensityComponentNames = getLogDensityComponentNames(object)
                    ,blockNames = getBlockNames(object)
                    ,nChain = nChainInPop
                    ,initialPopulation = object@initialPopulation[[iPop]]
                )
            })
        # TODO: check that enough samples (m0) in initial population
        object
    }
)

if(!exists("nSamplePopulations<-")) setGeneric("nSamplePopulations<-", function(object,value) standardGeneric("nSamplePopulations<-"))
setReplaceMethod("nSamplePopulations", signature=c("SampleLogs", "integer"), 
    function(object, value) {
        # adjust the sample space in population logs
        nSamplePopuplations <- if( length(value) == 1) rep(value, getNPopulation(object)) else value 
        if( length(nSamplePopuplations) != getNPopulation(object) ) stop("Need to provide number of samples for each population.")
        object@populationLogs <- mapply( function(populationLog, nSample){
                nSample(populationLog) <- nSample
                populationLog
            },object@populationLogs, nSamplePopuplations)
        object@sampleDimensions <- modifySampleDimensionsImpl(object@sampleDimensions
            ,nSamplePopulations = nSamplePopuplations)
        object
    })



if(!exists("getParametersForPopulation")) setGeneric("getParametersForPopulation", function(object, iPopulation) standardGeneric("getParametersForPopulation"))
#' @export
setMethod("getParametersForPopulation", signature="SampleLogs", 
    function( object, iPopulation){
        getParameters( object@populationLogs[[iPopulation]] )
    }
)

if(!exists("getParametersForAllChains")) setGeneric("getParametersForAllChains", function(object, iPopulations) standardGeneric("getParametersForAllChains"))
#' @export
setMethod("getParametersForAllChains", signature="SampleLogs", 
        function( object
            , iPopulations=1:getNPopulation(object)     ##<< the populations, for which to extract the parameters, by default all 
        ){
            if( missing(iPopulations) ) iPopulations=1:getNPopulation(object) 
            # first bind all the populations by last dim, i.e. chains
            unStacked <- do.call(abind, c(lapply( iPopulations, function(iPopulation){
                getParameters( object@populationLogs[[iPopulation]] )
            }), list(rev.along=1L)) )
            # next stack chains
            structure(unStacked, dim=c(nrow(unStacked), dim(unStacked)[2]*dim(unStacked)[3]), dimnames=dimnames(unStacked)[1:2])
        }
)


if(!exists("getLogDensityComponentsForPopulation")) setGeneric("getLogDensityComponentsForPopulation", function(object, iPopulation) standardGeneric("getLogDensityComponentsForPopulation"))
#' @export
setMethod("getLogDensityComponentsForPopulation", signature="SampleLogs", 
    function( object, iPopulation){
        getLogDensityComponents( object@populationLogs[[iPopulation]] )
    }
)

if(!exists("getParametersAndInitialForPopulation")) setGeneric("getParametersAndInitialForPopulation", function(object, iPopulation) standardGeneric("getParametersAndInitialForPopulation"))
#' @export
setMethod("getParametersAndInitialForPopulation", signature="SampleLogs", 
    function( object
        ### get the parameters including rows of initial values for given population index.
        , iPopulation   ##<< integer scalar: index of population
    ){
        ##value<< numeric matrix (nParm, nStepInitial+nStep, nChain)
        getParametersAndInitial( object@populationLogs[[iPopulation]] )
    }
)

if(!exists("setSampleRecord")) setGeneric("setSampleRecord", function(object, ... ) standardGeneric("setSampleRecord"))
setMethod("setSampleRecord", signature=c("SampleLogs"), 
    function( object
        , iSample
        , parameters                    ##<< numeric matrix (nParm,nPopulation*nChainInPopulation)
        , logDensityComponents          ##<< numeric matrix (nParm,nPopulation*nChainInPopulation)
        , stepTemperature=rep(1,getNPopulation(object)) ##<< numeric vector (nPopulation)
    ) {
        if( length(stepTemperature) != getNPopulation(object)) stop("Need to set temperature for each population.")
        for(iPopulation in 1:getNPopulation(object) ){
            iChains <- getIChainsForPopulation(object, iPopulation)
            object@populationLogs[[iPopulation]] <- setSampleRecord( object@populationLogs[[iPopulation]]
                , iSample=iSample
                ,parameters = parameters[,iChains ,drop=FALSE]
                ,logDensityComponents = logDensityComponents[,iChains ,drop=FALSE]
                ,stepTemperature = stepTemperature[iPopulation]
            )
        }
        object
    }
)

if(!exists("recordSampleLogChains")) setGeneric("recordSampleLogChains", function(object, sampleLogChains, ... ) standardGeneric("recordSampleLogChains"))
setMethod("recordSampleLogChains", signature=c("SampleLogs","list"), 
    function( object
        , sampleLogChains                       ##<< list of SampleLogChainImpl that take part in the step
        , iPopulations=1:getNPopulation(object) ##<< integer vector of population indices 
        , iSample0=rep(0L,getNPopulation(object))   ##<< integer vector (nPopulation): the step after which to record the sampleLog 
    ) {
        if( length(iSample0) != getNPopulation(object)) stop("need to specify iSample0 for each populations.")
        if( length(sampleLogChains) != length(iPopulations)*getNChainInPopulation(object))
            stop("number of chains in sampleLog does not match with argument iPopulations.")
        for( i in seq_along(iPopulations) ){
            iPop <- iPopulations[i]
            sampleLogPop <- object@populationLogs[[iPop]]
            # note that using i instead of iPop here, because iChain within chains in step
            iChainsInSamplesChains <- getIChainsForPopulation(object, i) 
            object@populationLogs[[iPop]] <- recordSampleLogs( object@populationLogs[[iPop]],
                sampleLogChains=sampleLogChains[iChainsInSamplesChains], iSample0=iSample0[iPop])
        }
        object
    }
)

if(!exists("getNInitialSamplePops")) setGeneric("getNInitialSamplePops", function(object) standardGeneric("getNInitialSamplePops"))
#' @export
setMethod("getNInitialSamplePops", signature="SampleLogs", function(object) {sapply(object@populationLogs, getNInitialSample)})

if(!exists("getInitialParameters")) setGeneric("getInitialParameters", function(object) standardGeneric("getInitialParameters"))
#' @export
setMethod("getInitialParameters", signature="SampleLogs", 
    function(object) {
        #populationLog <- object@populationLogs[[1]]
        do.call( cbind, sapply(object@populationLogs, function(populationLog){
                    adrop(getParametersAndInitial(populationLog)[,getNInitialSample(populationLog), ,drop=FALSE],2)
                }, simplify=FALSE)) 
    })

if(!exists("getSampleLogChainsForPopulations")) setGeneric("getSampleLogChainsForPopulations", function(object, iPopulations) standardGeneric("getSampleLogChainsForPopulations"))
setMethod("getSampleLogChainsForPopulations", signature=c("SampleLogs","integer"), 
    function(object
        , iPopulations
    ) {
        do.call(c, lapply(object@populationLogs[iPopulations], getSampleLogChains, blockNames=getBlockNames(object)) )            
    })



if(!exists("stepTemperaturesPopulations<-")) setGeneric("stepTemperaturesPopulations<-", function(object,...,value) standardGeneric("stepTemperaturesPopulations<-"))
setReplaceMethod("stepTemperaturesPopulations", signature=c("SampleLogs"), 
    function(object, iSample0, value) {
        if(length(value) != getNPopulation(object) ) stop("Mismatch in number of populations, when setting temperatures.")
        if( !length(object@populationLogs) )
            object <- initializeSampleLogs(object)
        object@populationLogs <- mapply("stepTemperatures<-", object@populationLogs, iSample0, value=value)                         
        object
    })
if(!exists("getStepTemperaturesPopulations")) setGeneric("getStepTemperaturesPopulations", function(object) standardGeneric("getStepTemperaturesPopulations"))
#' @export
setMethod("getStepTemperaturesPopulations", signature="SampleLogs", function(object) {  lapply(object@populationLogs, getStepTemperatures )})

if(!exists("getEndTemperaturesPopulations")) setGeneric("getEndTemperaturesPopulations", function(object) standardGeneric("getEndTemperaturesPopulations"))
#' @export
setMethod("getEndTemperaturesPopulations", signature="SampleLogs", 
    function(object) { 
        # populationLog <- object@populationLogs[[1]]
        do.call(cbind, lapply(object@populationLogs, function(populationLog){
                getStepTemperatures(populationLog)[,getNSample(populationLog)]
            }))
        ##value<< numeric matrix (nLogDensityComponent, nPopulation)
    })

if(!exists("setEndTemperaturesPopulations")) setGeneric("setEndTemperaturesPopulations", function(object, ...) standardGeneric("setEndTemperaturesPopulations"))
#' @export
setMethod("setEndTemperaturesPopulations", signature="SampleLogs", 
        function(object
            ,TEnd ##<< numeric matrix (nLogDensityComponent, nPopulation)
        ) {  
            if( !is.matrix(TEnd) ) stopDemac("setEndTemperaturesPopulations: Need to provide a matrix with TEnd.")
            if( nrow(TEnd) != getNLogDensityComponent(object)) stopDemac("setEndTemperaturesPopulations: number of rows in TEnd must match number of logDensityComponents.")
            if( ncol(TEnd) != getNPopulation(object)) stopDemac("setEndTemperaturesPopulations: nuber of columns in TEnd must match number of populations.")
            object@populationLogs <- lapply(seq_along(object@populationLogs), function(iPop){
                        setEndTemperature(object@populationLogs[[iPop]], TEnd[,iPop] )                   
                    })
            object
        })


if(!exists("stackChainsInPopulation")) setGeneric("stackChainsInPopulation", function(object,...) standardGeneric("stackChainsInPopulation"))
#' @export
setMethod("stackChainsInPopulation", signature="SampleLogs", 
    function(object
        , mergeMethod="stack"	##<< method mixing/merging the chains (see twMergeSequences)
    ){
        object@initialPopulation <- list()
        object@populationLogs <- lapply(object@populationLogs, stackChains, mergeMethod=mergeMethod)
        sDim <- object@sampleDimensions 
        object@sampleDimensions <- modifySampleDimensionsImpl(object@sampleDimensions
            ,nSamplePopulations = getNSamplePopulations(sDim)*getNChainInPopulation(sDim)
            ,nChainInPopulation = 1L
        )
        object
    })
if(!exists("asMcmc")) setGeneric("asMcmc", function(object,...) standardGeneric("asMcmc"))
#' @export
setMethod("asMcmc", signature="SampleLogs", 
        function(object
            , iParms=TRUE   ##<< index to parameter array (boolean, integer or names)
        ) {
            squeezedLogs <- if( getNPopulation(object) == 1L){ 
                        object 
                    } else {
                        squeeze(object, 1) # make same length
                    }  
            #nSample <- min(getNSamplePopulations(object) )
            mcmcs <- do.call(c, lapply(squeezedLogs@populationLogs, function(populationLog){
                                #parameters <- getParameters(populationLog)[,1:nSample, ,drop=FALSE]
                                #apply( parameters, 3, function(parametersChain){  list(mcmc(t(parametersChain))) })
                                apply( getParameters(populationLog), 3, function(parametersChain){  list(mcmc(t(parametersChain[iParms,]))) })
                            }))
            ##value<< mcmc.list of all chains, squeezed to match the population with shortest chains
            res <- mcmc.list(unlist(mcmcs, recursive=FALSE))
        })

if(!exists("squeeze")) setGeneric("squeeze", function(object,...) standardGeneric("squeeze"))
#' @export
setMethod("squeeze", signature="SampleLogs", 
        function(
                ### subsample population to nSample, defaulting to fraction*min(nSamplePop)  
                object
                , fraction=1.0  ##<< fraction of resulting samples per original number of samples 
                    ## (of population with smallest sample number)
                , nSample = round(min(getNSamplePopulations(object))*fraction)   ##<< integer vector (nPopulation): 
                    ## number of samples for each population
                    ## If lenght is 1, then this length is applied to all populations
        ) {
            ##details<< 
            ## By default all populations are subsampled, to match the sample number of shortest population
            nSampleOrig <- getNSamplePopulations(object)
            nSamplePopulations <- if( length(nSample)==1L ) rep(nSample, getNPopulation(object) ) else nSample
            if( length(nSamplePopulations) != getNPopulation(object) ) stop("Length of nSample must match the number of populations.")
            if( !all(nSample <= nSampleOrig)) stop("Must ensure that (nSample <= getNSamplePopulations(object))")
            iSamplesPops <- lapply( 1:getNPopulation(object), function(iPop){
                        round(seq(1L, nSampleOrig[iPop], length.out=nSamplePopulations[iPop]))
                    })
            subSample(object, iSamplesPops)
        })

if(!exists("subSample")) setGeneric("subSample", function(object,...) standardGeneric("subSample"))
#' @export
setMethod("subSample", signature="SampleLogs", 
        function(
                ##<< Retain the samples at given steps.
                object
                , iSamplesPops  ##<< list (nPopulation) of integer vector of sample indices.
                    ## If length is one, it is applied to all populations.
        ) {
            # one indices vector for each population
            if( !is.list(iSamplesPops) ) iSamplesPops <- rep( list(iSamplesPops), getNPopulation(object))
            if( length(iSamplesPops) != getNPopulation(object)) stop("Need to provide sample indices for each population.")
            object@initialPopulation <- list()
            object@populationLogs <- mapply(subSample, object@populationLogs, iSamplesPops, SIMPLIFY=FALSE)
            object@sampleDimensions <- modifySampleDimensionsImpl(object@sampleDimensions
                    ,nSamplePopulations = sapply(iSamplesPops, length)
            )
            object
        })


if(!exists("subsetTail")) setGeneric("subsetTail", function(object,...) standardGeneric("subsetTail"))
#' @export
setMethod("subsetTail", signature="SampleLogs", 
    function(
            ### Retain the \code{fraction} samples at the tail of the chains
            object      
            , fraction  ##<< numeric scalar: the fraction of samples to retain
            , nSample = round(getNSamplePopulations(object)*fraction)   ##<< alternative way of specifying the number of samples in the tails
    ) {
        if( !missing(fraction) && (fraction < 0 || fraction > 1)) stopDemac("Proportion must be within 0 and 1.")
        nSample0 <- getNSamplePopulations(object)
        nPop <- length(nSample0)
        if( length(nSample)==1L) nSample <- rep(nSample, nPop)
        if( length(nSample) != nPop) stopDemac("length of nSample must correspond to number of popultions.")
        iShort <- which(nSample0 < nSample)
        if( length(iShort) ){
            warning("SampleLogs.subsetTail: requesting more samples than available, returning all available samples.")
            nSample[iShort] <- nSample0[iShort]                
        } 
        iSamplesPops <- lapply( 1:nPop, function(iPop){
                (nSample0[iPop]+1L-nSample[iPop]):nSample0[iPop]
            })
        subSample(object, iSamplesPops)
    })


if(!exists("subsetPopulations")) setGeneric("subsetPopulations", function(object,...) standardGeneric("subsetPopulations"))
#' @export
setMethod("subsetPopulations", signature="SampleLogs", 
        function(
                ### Retain the \code{fraction} samples at the tail of the chains
                object      
                , iPopulations=1L  ##<< index (boolean, integer) to populations 
        ) {
            # one indices vector for each population
            object@sampleDimensions <- modifySampleDimensionsImpl(object@sampleDimensions
                    ,nSamplePopulations = getNSamplePopulations(object@sampleDimensions)[iPopulations]
            )
            object@initialPopulation <- object@initialPopulation[iPopulations]
            object@populationLogs <- object@populationLogs[iPopulations]
            object
        })



if(!exists(".appendSampleLogs")) setGeneric(".appendSampleLogs", function(object,...) standardGeneric(".appendSampleLogs"))
setMethod(".appendSampleLogs", signature="SampleLogs", 
    function(object, sampleLogPopulations) {
        if( length(sampleLogPopulations) != getNPopulation(object) ) stop("Need to provide a sampleLog for each population when appending.")
        object@populationLogs <- mapply( function(log1,log2){ appendSampleLog(log1,log2)}
            , object@populationLogs, sampleLogPopulations)
        object@sampleDimensions <- modifySampleDimensionsImpl(object@sampleDimensions
            , nSamplePopulations = getNSamplePopulations(object@sampleDimensions) + sapply(sampleLogPopulations,getNSample) )
        object    
    })

if(!exists("applyPopulations")) setGeneric("applyPopulations", function(object,...) standardGeneric("applyPopulations"))
#' @export
setMethod("applyPopulations", signature="SampleLogs", 
    function(
            ### lapply a function to all popultions
        object  ##<< SampleLogs object
        ,FUN    ##<< function to be applied. 
            ##<< Signature must correspond to FUN (SampleLogs, indexOfPopultion, ... )
            ##<< e.g. \code{\link{getParametersForPopulation}}
        ,...    ##<< further arguments to FUN
    ) {
        lapply(1:getNPopulation(object),function(iPop){
                FUN(object, iPop, ...)
            })
        ##value<< a list with results of applying FUN to each population
    })

if(!exists("computeTemperatedLogDensityForPopulation")) setGeneric("computeTemperatedLogDensityForPopulation", function(object,...) standardGeneric("computeTemperatedLogDensityForPopulation"))
#' @export
setMethod("computeTemperatedLogDensityForPopulation", signature="SampleLogs", 
    function(object
        ### multiply logDensityComponents by given temperatures
        , iPopulation                               ##<< scalar integer: population
        ,temperature = getEndTemperaturesPopulations(object)[,iPopulation]  ##<< numeric vector (nResComp): temperatures to multiply logDensityComponents
            ## Also a matrix specifying a temperature for each step is possible (\code{\link{getStepTemperatures(object)}})
    ) {
        computeTemperatedLogDensity( object@populationLogs[[iPopulation]], temperature=temperature)
    })

if(!exists("computeSampleRanksForPopulation")) setGeneric("computeSampleRanksForPopulation", function(object,...) standardGeneric("computeSampleRanksForPopulation"))
#' @export
setMethod("computeSampleRanksForPopulation", signature="SampleLogs", 
    function(object
        ### rank the samples 
        , iPopulation   ##<< scalar integer: index of the population of interest
    ) {
        ##value<<
        computeSampleRanks( object@populationLogs[[iPopulation]], object)
    })

if(!exists("blockDimensions<-")) setGeneric("blockDimensions<-", function(object,value) standardGeneric("blockDimensions<-"))
setReplaceMethod("blockDimensions", signature=c("SampleLogs", "BlockDimensions"), function(object,value){ 
            blockDimensions(object@sampleDimensions) <- value
            object
        }) 

if(!exists("getLogDensitiesForPopulation")) setGeneric("getLogDensitiesForPopulation", function(object,...) standardGeneric("getLogDensitiesForPopulation"))
#' @export
setMethod("getLogDensitiesForPopulation", signature="SampleLogs", 
        function(object
            ### Get the sum of temperated logDensityComponents within blocks
            , iPopulation   ##<< scalar integer: population 
            ,temperature = getEndTemperaturesPopulations(object)[,iPopulation]  ##<< numeric vecotr (nResComp) of temperatures, that logDensityComponents are multiplied before sum across blocks
        ) {
            temperatedLogDensityComponents <- computeTemperatedLogDensityForPopulation(object, iPopulation, temperature=temperature)
            logDenBlocksChains <- lapply( 1:getNChain(object), function(iChain){
                        .sumLogDensitiesAcrossBlocks(adrop(temperatedLogDensityComponents[,,iChain ,drop=FALSE],3), object@sampleDimensions) 
                    })
            ##value<< numeric array (nBlock, nSample, nChain) of block logDensity  
            abind(logDenBlocksChains, rev.along=0)
        })

#------- delegation of SampleDimensions interface to SampleDimensionsImpl
#' @export
setMethod("getNBlock", signature="SampleLogs", function(object) {getNBlock(object@sampleDimensions)})
#' @export
setMethod("getBlockNames", signature="SampleLogs", function(object) {getBlockNames(object@sampleDimensions)})

#' @export
setMethod("getNParameter", signature="SampleLogs", function(object) {getNParameter(object@sampleDimensions)})
#' @export
setMethod("getNParameterWithProposal", signature="SampleLogs", function(object) {getNParameterWithProposal(object@sampleDimensions)})
#' @export
setMethod("getParameterNames", signature="SampleLogs", function(object) {getParameterNames(object@sampleDimensions)})
#' @export
setMethod("getIParametersUpdatedByBlock", signature=c("SampleLogs","integer"), function(object, blockIndex) {getIParametersUpdatedByBlock(object@sampleDimensions, blockIndex)})
#' @export
setMethod("getIParametersUpdatedByBlock", signature=c("SampleLogs","character"), function(object, blockIndex) {getIParametersUpdatedByBlock(object@sampleDimensions, blockIndex)})

#' @export
setMethod("getNLogDensityComponent", signature="SampleLogs", function(object) {getNLogDensityComponent(object@sampleDimensions)})
#' @export
setMethod("getLogDensityComponentNames", signature="SampleLogs", function(object) {getLogDensityComponentNames(object@sampleDimensions)})
#' @export
setMethod("getILogDensityComponentsByBlock", signature=c("SampleLogs","integer"), function(object, blockIndex) {getILogDensityComponentsByBlock(object@sampleDimensions, blockIndex)})
#' @export
setMethod("getILogDensityComponentsByBlock", signature=c("SampleLogs","character"), function(object, blockIndex) {getILogDensityComponentsByBlock(object@sampleDimensions, blockIndex)})

#' @export
setMethod("getNSamplePopulations", signature="SampleLogs", function(object) {getNSamplePopulations(object@sampleDimensions)})
#' @export
setMethod("getNPopulation", signature="SampleLogs", function(object) {getNPopulation(object@sampleDimensions)})
#' @export
setMethod("getNChain", signature="SampleLogs", function(object) {getNChain(object@sampleDimensions)})
#' @export
setMethod("getNChainInPopulation", signature="SampleLogs", function(object) {getNChainInPopulation(object@sampleDimensions)})
# ' @export
# setMethod("getNChainInChainGroup", signature="SampleLogs", function(object) {getNChainInChainGroup(object@sampleDimensions)})
#' @export
setMethod("getNChainGroupsInPopulation", signature="SampleLogs", function(object) {getNChainGroupsInPopulation(object@sampleDimensions)})
#' @export
setMethod("getIChainsForPopulation", signature=c("SampleLogs","integer"), function(object, iPopulation) {getIChainsForPopulation(object@sampleDimensions, iPopulation)})

#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("SampleLogs")
if(!exists("getSampleDimensions")) setGeneric("getSampleDimensions", function(object) standardGeneric("getSampleDimensions"))
#' @export
setMethod("getSampleDimensions", signature="SampleLogs", function(object) {object@sampleDimensions})

if(!exists("getInitialPopulation")) setGeneric("getInitialPopulation", function(object) standardGeneric("getInitialPopulation"))
#' @export
setMethod("getInitialPopulation", signature="SampleLogs", function(object) {object@initialPopulation})

if(!exists("getSampleLogOfPopulation")) setGeneric("getSampleLogOfPopulation", function(object,...) standardGeneric("getSampleLogOfPopulation"))
#' @export
setMethod("getSampleLogOfPopulation", signature="SampleLogs", function(object, iPopulation) {object@populationLogs[[iPopulation]]})


