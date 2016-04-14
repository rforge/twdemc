#' @include SampleLogChainImpl.R

#' @export
setClass("SampleLogPopulation",
    representation(  
        nInitialSample = "integer"	           ##<< scalar integer: number of samples in initial population
        ,parameters = "array"                  ##<< numeric array (nParm, nInterval, nChain): sampled parameters
        ,logDensityComponents = "array"        ##<< numeric array (nResComp, nInterval, nChain): result components of logDensity for those parameters
        ,stepTemperatures = "matrix"	       ##<< numeric matrix (nResComp, nInterval): temperatures at given step (all chains are at the same temperature)
    ),
    prototype=list(nInitialSample=0L)
)

setMethod("show", "SampleLogPopulation",
    function(object){
        cat("(nSample=",getNSample(object)
            ," nChain=",dim(object@parameters)[3]
            ," nParm=",nrow(object@parameters)
            ," nLogDensityComponent=",nrow(object@logDensityComponents)
            ,sep="")
        cat(")\n")
    })

#TODO record from chains, with keeping track on nSample

if(!exists("initializeLog")) setGeneric("initializeLog", function(object,...) standardGeneric("initializeLog"))
setMethod("initializeLog", signature="SampleLogPopulation", 
    function(object,
        ### allocate and initialize storage (parameters, logDensityComponents, countAcceptedInterval) 
        nSample,                    ##<< scalar integer: number of samples, i.e. rows, to store
        parameterNames,              
        logDensityComponentNames,   
        blockNames,
        nChain,
        initialPopulation=numeric(0)
    ) {
        #if( missing(initialPopulation) ) initialPopulation=numeric(0)
        object@nInitialSample <- if(length(initialPopulation)) ncol(initialPopulation) else 0L
        object@parameters <- array(NA_real_, dim=c(length(parameterNames),object@nInitialSample+nSample,nChain), 
            dimnames=list(parameterNames,NULL,NULL))
        if( object@nInitialSample != 0L)
            object@parameters[,1:object@nInitialSample,] <- initialPopulation    
        object@logDensityComponents <- array(NA_real_, dim=c(length(logDensityComponentNames),nSample,nChain), 
            dimnames=list(logDensityComponentNames,NULL,NULL))
        object@stepTemperatures <- matrix(1L, length(logDensityComponentNames), ncol=nSample,
            dimnames=list(logDensityComponentNames,NULL))
        object
    }
)

if(!exists("initialPopulation<-")) setGeneric("initialPopulation<-", function(object,value) standardGeneric("initialPopulation<-"))
setReplaceMethod("initialPopulation", signature=c("SampleLogPopulation", "array"), 
    function(object
        ,value
    ){
        parameterNames <- rownames(object@parameters)
        object@parameters <- abind(value, getParameters(object), along=2 )
        object@nInitialSample <- if(length(value)) ncol(value) else 0L
        rownames(object@parameters) <- parameterNames
        object
    })

if(!exists("getNSample")) setGeneric("getNSample", function(object) standardGeneric("getNSample"))
#' @export
setMethod("getNSample", signature="SampleLogPopulation", function(object) { ncol(object@logDensityComponents) })

if(!exists("getNChain")) setGeneric("getNChain", function(object) standardGeneric("getNChain"))
#' @export
setMethod("getNChain", signature="SampleLogPopulation", function(object) { dim(object@parameters)[3]} )

if(!exists("getSampleLogChains")) setGeneric("getSampleLogChains", function(object, ...) standardGeneric("getSampleLogChains"))
setMethod("getSampleLogChains", signature="SampleLogPopulation", function(object, blockNames) {
        lapply(1:getNChain(object), function(iChain){
                new("SampleLogChainImpl", nSample=getNSample(object)
                    , parameters=adrop(getParameters(object)[,,iChain ,drop=FALSE],3)
                    , logDensityComponents=adrop(object@logDensityComponents[,,iChain ,drop=FALSE],3)
                    , proportionAcceptedInInterval=matrix(0.25, nrow=length(blockNames), ncol=ncol(object@logDensityComponents),dimnames=list(blockNames,NULL))
                )
            })  
        
    })

if(!exists("recordSampleLogs")) setGeneric("recordSampleLogs", function(object, sampleLogChains, ... ) standardGeneric("recordSampleLogs"))
setMethod("recordSampleLogs", signature=c("SampleLogPopulation","list"), 
    function( object, sampleLogChains, iSample0 ){
        for( iChain in 1:getNChain(object) ){
            sampleLog <- sampleLogChains[[iChain]]
            nSample <- getNSample(sampleLog)
            object@parameters[,object@nInitialSample+iSample0+(1:nSample),iChain] <- getParameters(sampleLog)
            object@logDensityComponents[,iSample0+(1:nSample),iChain] <- getLogDensityComponents(sampleLog)
        }
        object
    }
)

if(!exists("setSampleRecord")) setGeneric("setSampleRecord", function(object, ... ) standardGeneric("setSampleRecord"))
setMethod("setSampleRecord", signature=c("SampleLogPopulation"), 
    function( object
        , iSample                       ##<< integer scalar: position at which to record       
        , parameters                    ##<< numeric matrix (nParm, nChain): parameters
        , logDensityComponents          ##<< numeric matrix (nLogDenComp, nChain): 
        , stepTemperature=1                 ##<< numeric scalar: current temperature 
    ){
        if( nrow(parameters) != nrow(object@parameters) ) stop("mismatch in number of parameters")
        object@parameters[,object@nInitialSample+iSample,] <- parameters
        object@logDensityComponents[,iSample,] <- logDensityComponents
        object@stepTemperatures[,iSample] <- stepTemperature
        object
    })

if(!exists("stepTemperatures<-")) setGeneric("stepTemperatures<-", function(object,..., value) standardGeneric("stepTemperatures<-"))
setReplaceMethod("stepTemperatures", signature=c("SampleLogPopulation"), 
    function(object
        ### Set the temperatures for each step
        , iSample0  ##<< integer scalar: index of record before the records for which to set temperatures 
        , value     ##<< numeric vector (nLogDensityComponents, nStep): temperature for each step
    ) {
        if( nrow(value) != nrow(object@logDensityComponents) ) stop("Mismatch in numbers of logDensityComponents, when setting temperatures.")
        if( ncol(value) != getNSample(object)-iSample0 ) stop("Mismatch in numbers of samples, when setting temperatures.")
        object@stepTemperatures[,iSample0+(1:ncol(value))] <- value
        object
    })

if(!exists("getEndTemperature")) setGeneric("getEndTemperature", function(object, ... ) standardGeneric("getEndTemperature"))
#' @export
setMethod("getEndTemperature", signature=c("SampleLogPopulation"), 
        function(object
                ### Get the temperature at the last step
        ) {
            ##value<< numeric vector (nResComp) of temperatures for each logDensityComponent at the end of the sampling
            object@stepTemperatures[,ncol(object@stepTemperatures)] 
        })

if(!exists("stackChains")) setGeneric("stackChains", function(object,...) standardGeneric("stackChains"))
#' @export
setMethod("stackChains", signature="SampleLogPopulation", 
    function(object
        ,mergeMethod="stack"	##<< method of mixing the several chains (see twMergeSequences)
    ) {
        object@parameters <- getParameters(object)  # drop initial population
        object@nInitialSample <- 0L
        nChain <- dim(object@parameters)[3]
        nSample <- dim(object@parameters)[2]
        parameters <- abind( lapply( 1:nChain, function(iChain){ object@parameters[,,iChain ,drop=FALSE]}), along=2 )
        logDensityComponents <- abind( lapply( 1:nChain, function(iChain){ object@logDensityComponents[,,iChain ,drop=FALSE]}), along=2 )
        stepTemperatures <- abind( lapply( 1:nChain, function(iChain){object@stepTemperatures}), along=2 )
        if( mergeMethod != "stack"){
            chainInd <- twMergeSequences( rep(nSample,nChain), mergeMethod=mergeMethod )
            for(iChain in 1:nChain){
                iPos <- which(chainInd==iChain)
                parameters[,iPos,1] <- object@parameters[,,iChain]
                logDensityComponents[,iPos,1] <- object@logDensityComponents[,,iChain]
                stepTemperatures[,iPos] <- object@stepTemperatures 
            }
        }
        object@parameters <- parameters
        object@logDensityComponents <- logDensityComponents
        object@stepTemperatures <- stepTemperatures
        object
    })

if(!exists("subSample")) setGeneric("subSample", function(object,...) standardGeneric("subSample"))
#' @export
setMethod("subSample", signature="SampleLogPopulation", 
    function(object
        ,iSamples
    ) {
        if( missing(iSamples) ) stop("Need to provide argument iSample.")
        if( is.numeric(iSamples) && max(iSamples) > getNSample(object)) stop("tried to subset non-existing sample indices.")
        object@parameters <- getParameters(object)  # drop initial population
        object@nInitialSample <- 0L
        #
        object@parameters <- object@parameters[,iSamples, ,drop=FALSE]
        object@logDensityComponents <- object@logDensityComponents[,iSamples, ,drop=FALSE]
        object@stepTemperatures <- object@stepTemperatures[,iSamples ,drop=FALSE]
        object
    })

if(!exists("appendSampleLog")) setGeneric("appendSampleLog", function(object,...) standardGeneric("appendSampleLog"))
setMethod("appendSampleLog", signature="SampleLogPopulation", 
    function(object
        ,sampleLog
    ) {
        if( getNChain(sampleLog) != getNChain(object)) stop("mismatch in chain number when appending sampleLog.")
        #
        object@parameters <- abind( getParameters(object), getParameters(sampleLog), along=2)
        object@logDensityComponents <- abind( getLogDensityComponents(object), getLogDensityComponents(sampleLog), along=2)
        object@stepTemperatures <- abind( getStepTemperatures(object), getStepTemperatures(sampleLog), along=2)
        object@nInitialSample <- 0L
        object
    })

if(!exists("nSample<-")) setGeneric("nSample<-", function(object,value) standardGeneric("nSample<-"))
setReplaceMethod("nSample", signature=c("SampleLogPopulation", "integer"), 
    function(object, value) {
        if( value < getNSample(object) ){
            # when nSample is smaller than actual, subset
            object@parameters <- object@parameters[,1:(object@nInitialSample+value), ,drop=FALSE] 
            object@logDensityComponents <- object@logDensityComponents[,1:value, ,drop=FALSE]
            object@stepTemperatures <- object@stepTemperatures[,1:value ,drop=FALSE]
        } else if( value > getNSample(object) ){
            # when new nSample is larger than current, append empty array
            dnSample <- value - getNSample(object)
            nChain <- getNChain(object)
            object@parameters <- abind( object@parameters
                , array(NA_real_, dim=c(nrow(object@parameters),dnSample,nChain))
                , along=2)
            object@logDensityComponents <- abind( object@logDensityComponents
                , array(NA_real_, dim=c( nrow(object@logDensityComponents),dnSample,nChain))
                , along=2)
            object@stepTemperatures <- abind( object@stepTemperatures
                , matrix(1L, nrow=nrow(object@stepTemperatures), ncol=dnSample)
                , along=2)
        }
        object
    })

if(!exists("computeTemperatedLogDensity")) setGeneric("computeTemperatedLogDensity", function(object,...) standardGeneric("computeTemperatedLogDensity"))
#' @export
setMethod("computeTemperatedLogDensity", signature="SampleLogPopulation", 
    function(object
        ### devide logDensityComponents by given temperature
        ,temperature = getEndTemperature(object)  ##<< numeric matrix (nResComp, nStep): see \code{\link{getStepTemperatures(object)}}.
            ## If a vector is provided, then it is used for all steps
    ) {
        logDen <- getLogDensityComponents(object)
        if( !is.matrix(temperature)){
            logDen / temperature  # T is recycled in multiplication
        } else {
            # MAYBE: check if abind is needed or if usual recycling is sufficient
            T <- abind( rep(list(temperature), getNChain(object)), rev.along=0)
            logDen / T
        }  
    })


#' @export
computeBlockLogDensities <- function(
        ### sum temperated logDensities for blocks
        object                  ##<< SampleLogPopulation object
        ,blockDimensions        ##<< object specifying which densities are used by which blocks
        ,temperature = getEndTemperature(object)  ##<< numeric matrix (nResComp, nStep) or numeric vector (nResComp). (see computeTemperatedLogDensity)
){
    logDenT <- computeTemperatedLogDensity(object, temperature=temperature)
    if( nrow(logDenT)==0L ){
        #recover()
        warning("computeSampleRanks: no logDensities provided.")  
    } 
    chainBlockLogDensities <- lapply( 1:getNChain(object), function(iChain){
                logDenTc <- adrop(logDenT[,,iChain ,drop=FALSE],3)
                logDenBlocks <- .sumLogDensitiesAcrossBlocks(logDenTc, blockDimensions)
            })
    ##value<< numeric array (nBlock x nStep x nChain)
    ans <- abind( chainBlockLogDensities, rev.along=0 )
    ans    
}


if(!exists("computeSampleRanks")) setGeneric("computeSampleRanks", function(object,...) standardGeneric("computeSampleRanks"))
#' @export
setMethod("computeSampleRanks", signature="SampleLogPopulation", 
    function(object
        , blockDimensions
        ,temperature = getEndTemperature(object)  ##<< numeric matrix (nResComp, nStep) or numeric vector (nResComp).
        ,decreasing=TRUE          ##<< if TRUE start with the highest ranks, i.e best samples 
    ) {
        logDenBlocks <- computeBlockLogDensities(object, blockDimensions=blockDimensions, temperature=temperature)
        #aaply fails on providing logDenBlocks with a degenerate dimesion in 1L or 2L
        chainRanks <- apply( logDenBlocks, 3L, function(logDenBlocksChain){
                    .rankBlockLogDensities( logDenBlocksChain, decreasing = decreasing )
                })
        ##value<< integer matrix (nSample x nChain) specifying the sample number of rank at vector position
        chainRanks
    })

.sumLogDensitiesAcrossBlocks <-  function(
    ### sum logDensity components across blocks
    logDensityComponents	##<< numeric matrix (nLogensitiyComponents x nStep)
    ,blockDimensions    ##<< Object of class BlockDimensions
){
    ##seealso<<
    ## \code{\link{orderLogDen}}, \code{\link{calcTemperatedLogDen}}
    iPosBlocks <-  lapply( 1:getNBlock(blockDimensions), function(iBlock){ getILogDensityComponentsByBlock(blockDimensions,iBlock) })
    iBlocksWithLogDen <- which( sapply(iPosBlocks, length ) != 0 ) 
    logDensityComponentsBlock <-
        if( length(iBlocksWithLogDen) == 0){
            logDensityComponentsBlock <- logDensityComponents # should be already zero row matrix 
        } else if( length(iBlocksWithLogDen) == 1){
            logDensityComponentsBlock <- matrix(colSums(logDensityComponents), nrow=1
                , dimnames=list(getBlockNames(blockDimensions)[iBlocksWithLogDen], NULL))
        } else {
            #iBlock = 1
            logDensityComponentsBlock <- do.call( rbind, tmp <- lapply( iBlocksWithLogDen, function(iBlock){
                        logDen <- colSums(logDensityComponents[iPosBlocks[[iBlock]], ,drop=FALSE])
                    }))
            rownames(logDensityComponentsBlock) <- getBlockNames(blockDimensions)[iBlocksWithLogDen]
            logDensityComponentsBlock
        }
    ##value<< numeric matrix (nBlock, nStep): of sum across logDensities
}


.rankBlockLogDensities <- function(
    ### rank the log-Densities, starting with highest densities
    logDensityComponentsBlock          ##<< numeric matrix ( nBlockWithLogDen, nStep ) 
    ,decreasing=TRUE          ##<< argument to order, to start with the highest ranks, i.e best models
){
    nBlock <- nrow(logDensityComponentsBlock)
    ##seealso<<
    ## \code{\link{computeSampleRanks}}
    ## \code{\link{.sumLogDensitiesAcrossBlocks}}
    #        
    sampleRanks <- if( nBlock==0) {
            # when no logDensities are given, return the current order
            1:ncol(logDensityComponentsBlock)
        } else if( nBlock==1){
            order( logDensityComponentsBlock, decreasing=decreasing )  
        }else{
            # calculate ranks for each logDensityComponent
            r <- apply( logDensityComponentsBlock, 1, rank )
            # get the minimum rank across components
            minR <- apply(r,1,min)
            # sort by this minimum rank
            order(minR, decreasing=decreasing)
        }
    ##value<< integer vector of indices (, ordered by ranks of the samples
    sampleRanks
}

# implementation of SampleLog

if(!exists("getParametersAndInitial")) setGeneric("getParametersAndInitial", function(object) standardGeneric("getParametersAndInitial"))
setMethod("getParametersAndInitial", signature="SampleLogPopulation", function(object) {object@parameters})


if(!exists("getParameters")) setGeneric("getParameters", function(object) standardGeneric("getParameters"))
#' @export
setMethod("getParameters", signature="SampleLogPopulation", function(object) {object@parameters[,(object@nInitialSample+1):ncol(object@parameters),, drop=FALSE]})
if(!exists("getLogDensityComponents")) setGeneric("getLogDensityComponents", function(object) standardGeneric("getLogDensityComponents"))
#' @export
setMethod("getLogDensityComponents", signature="SampleLogPopulation", function(object) {
            object@logDensityComponents
        })

if(!exists("getParametersForSampleIndex")) setGeneric("getParametersForSampleIndex", function(object,iSample) standardGeneric("getParametersForSampleIndex"))
setMethod("getParametersForSampleIndex", signature="SampleLogPopulation", function(object,iSample) {
        adrop(object@parameters[,object@nInitialSample+iSample, ,drop=FALSE],2)})
if(!exists("getLogDensityComponentsForSampleIndex")) setGeneric("getLogDensityComponentsForSampleIndex", function(object, iSample) standardGeneric("getLogDensityComponentsForSampleIndex"))
setMethod("getLogDensityComponentsForSampleIndex", signature="SampleLogPopulation", function(object, iSample) {
        adrop(object@logDensityComponents[,iSample, ,drop=FALSE],2)})

if(!exists("getBestSamples")) setGeneric("getBestSamples", function(object, ...) standardGeneric("getBestSamples"))
#' @export
setMethod("getBestSamples", signature=c("SampleLogPopulation"),
        function( object
            ##<< get a singleChained log of best samples
                 ,blockDimensions    ##<< BlockDimensions object with information on blocks
                ,fractionBest=0.10   ##<< numeric scalar: proportion of the sample, defaults to best 10%
                ,nBest=round(getNSample(sampleLogStacked)*fractionBest) ##<< alternative way of specifying the numbe of best samples
        ) {
            sampleLogStacked <- stackChains(object) # combine all the chains to one big
            iBest <- computeSampleRanks(sampleLogStacked, blockDimensions)[1:nBest,1L]
            subSample(sampleLogStacked, iBest)
        })

if(!exists("setEndTemperature")) setGeneric("setEndTemperature", function(object, value) standardGeneric("setEndTemperature"))
setMethod("setEndTemperature", signature="SampleLogPopulation", function(object
    , value ##<< numeric vector (nLogDenComp): temperature of last sample for all logDensity components
) {
            object@stepTemperatures[,getNSample(object)] <- value
            object
        })



#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("SampleLogPopulation")
if(!exists("getNInitialSample")) setGeneric("getNInitialSample", function(object) standardGeneric("getNInitialSample"))
#' @export
setMethod("getNInitialSample", signature="SampleLogPopulation", function(object) {object@nInitialSample})

if(!exists("getStepTemperatures")) setGeneric("getStepTemperatures", function(object) standardGeneric("getStepTemperatures"))
#' @export
setMethod("getStepTemperatures", signature="SampleLogPopulation", function(object) {object@stepTemperatures})



