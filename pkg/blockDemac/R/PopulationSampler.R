#' @include SampleLogs.R
#' @include ChainSamplerImpl.R
#' @include JumpProposer.R
#' @include DEJumpProposer.R
#' @include AcceptanceTracker.R


#' @export
newPopulationSampler <- function(
        ### create a new PopulationSampler based on given blockSpecifications
        blockSpecifications ##<< list of blockSpecifications (see \code{\link{blockSpec}})
        ,thetaPrior		##<< numeric vector (nParm) of point estimate of the mean
        ,covarTheta 	##<< numeric matrix (nParm x nParm) the covariance of parameters.
        ,...            ##<< other arguments to \code{\link{createSampleLogsWithInitialPopulation}}, e.g. \code{Zinit}, \code{nChainGroupsInPopulation}
        ## such as \code{nPopulation} and \code{nChainInPopulation}
        ,intermediateSpecifications=list()       ##<< list of intermediateSpecifications (see \code{\link{intermediateSpec}})
        ,chainSamplerTemplate=new("ChainSamplerImpl")   ##<< subClass of ChainSampler that will be initialized and used
        ,className="PopulationSampler"
        ,subSpace=new("SubSpace")                       ##<< constraint of sampling range
        ,sampleLogs     ##<< alternative to specify initial population, instead of thetaPrior and covarTheta
){
    namesTheta <- if( !missing(sampleLogs) ){ getParameterNames(sampleLogs) } else
                  if( hasArg("Zinit") ){ rownames(list(...)$Zinit) } else
                  if( !missing(thetaPrior) ){ names(thetaPrior) } else
                  stopDemac("need specify initial population either by Zinit, sampleLogs, or (thetaPrior, covarTheta)")
    blockUpdaters <- newBlockUpdaters(blockSpecifications, intermediateSpecifications, namesTheta)
    chainSampler <- initializeBlockUpdaters(chainSamplerTemplate, blockUpdaters)
    sampler <- new(className, chainSampler=chainSampler)
    sampleLogs(sampler) <-if( missing(sampleLogs) ){
        initialLogs <- createSampleLogsWithInitialPopulation(sampler
                ,thetaPrior=thetaPrior
                ,covarTheta=covarTheta
                ,subSpace=subSpace
                ,...
        ) 
        initialLogs
    } else {
        sampleLogs
    }
    subSpacePopulations(sampler) <- rep(list(subSpace), getNPopulation(getSampleLogs(sampler))) 
    sampler            
}


setOldClass("cluster")
setClass("SOCKcluster", contains = "cluster")
setOldClass("SOCKcluster", S4Class = "SOCKcluster")

#' @export
setClass("PopulationSampler",
        representation(
                chainStates = "list"
                ,sampleLogs = "SampleLogs"
                ,chainSampler = "ChainSampler"
                ,jumpProposer = "JumpProposer"
                ,acceptanceTracker = "AcceptanceTracker"
                ,progressOutput = "character"
                ,subSpacePopulations = "list"	                    
                ,stepTemperaturesPopulations = "list"	 
                ,cluster = "cluster"
                ,nIntervalPerRange = "integer"	        ##<< number of thinning intervals before regenerating jumps
                ,nSampleBeforeSampling = "integer"	    # integer vector (nPopulation): number of samples in sampleLog, before sampling 
                ,nSampleBatchPopulations = "integer"	            # integer vector (nPopulation): number of samples sample and append to sampleLog
                ,isBurnin = "logical"	        ##<< if current sampling proceeds in burnin mode
                ,nRepresentativeRows= "integer"	        ##<< number of independent rows necessary to represent parameter space with current number of chainsInPopulation
        )
        ,prototype( 
                chainSampler = new("ChainSamplerImpl")
                , jumpProposer = new("DEJumpProposer")  
                , acceptanceTracker = new("AcceptanceTracker")
                , sampleLogs = new("SampleLogs")
                , cluster=structure(list(),class=c("SOCKcluster","cluster"))
                , progressOutput="."
                , nIntervalPerRange=2L
                , isBurnin=TRUE
                , nRepresentativeRows=NA_integer_
        )
)


if(!exists("createSampleLogsWithInitialPopulation")) setGeneric("createSampleLogsWithInitialPopulation", function(object,...) standardGeneric("createSampleLogsWithInitialPopulation"))
setMethod("createSampleLogsWithInitialPopulation", signature="PopulationSampler", 
        function( object
                ### Generate an initial population of states for \code{\link{twDEMCBlockInt}}.
                ,thetaPrior		##<< numeric vector (nParm) of point estimate of the mean
                ,covarTheta 	##<< numeric matrix (nParm x nParm) the covariance of parameters.
                ##<< Alternatively, can be given as vector (nParm) of the diagonal, i.e. variances, if all dimensions are independent
                , nPopulation=2L
                , nChainInPopulation=3L
                , nChainGroupsInPopulation=1L
                ,m0= computeNRepresentativeRows(length(thetaPrior),nChainInPopulation)/(m0FiniteFac)
                ### number of initial states for each chain
                ,m0FiniteFac=1	##<< use a factor smaller than 1 to increase default m0 to account for only a portion of proposal results in finite densities
                ,doIncludePrior=TRUE
                ### If TRUE, then set last sample of chain 1 to the prior estimate, 
                ### which might be already a kind of best estimates by an optimization.
                ,subSpace=new("SubSpace")                       ##<< constraint of sampling range
                ,Zinit          ##<< sample (nParm, nStep) as an alternative to (thetaPrior,covarTheta) to specify initial logs
        ){
            ##seealso<<  
            ## \code{\link{computeNRepresentativeRows}}
            #?rmvnorm #in package mvtnorm
            if( !missing(Zinit) ) thetaPrior <- Zinit[,1L]
            blockDimensions <- getBlockDimensions(object@chainSampler)
            if( getNBlock(blockDimensions)==0 ) stop("Need to setup chainSampler with blockSpecificiations before computing setting up sampleLog.")
            sampleDimensions <- initializeSampleDimensionsImpl(
                    new("SampleDimensionsImpl")
                    , blockDimensions = blockDimensions
                    , nSamplePopulations=rep(1L, nPopulation)
                    , nChainInPopulation=nChainInPopulation
                    , nChainGroupsInPopulation = nChainGroupsInPopulation
            )
            nChain <- getNChain(sampleDimensions) 	##<< number of chains to run
            nParameter <- getNParameter(sampleDimensions)
            parameterNames <- getParameterNames(sampleDimensions)
            if( length(thetaPrior) != nParameter )
                stop("Mismatch between length of prior and number of parameters in Sampler.")
            if( !length(names(thetaPrior)) ) names(thetaPrior) <- parameterNames 
            if( !all(names(thetaPrior) == parameterNames) )
                stop("Mismatch between parameternames in prior and in Sampler.")
            ZinitChains <- if( missing(Zinit) ) computeInitialPopulation(
                    thetaPrior, covarTheta
                    #, sampleDimensions=sampleDimensions
                    ,nChainInPopulation = getNChainInPopulation(sampleDimensions)
                    ,nPopulation = getNPopulation(sampleDimensions)
                    , m0=m0, m0FiniteFac=m0FiniteFac, doIncludePrior=doIncludePrior
                    , subSpace=subSpace
            ) else {
                # condense Zinit to m0 * nChain
                if( ncol(Zinit)  < (m0*nChain/2) ) stopDemac("too few initial samples")
                structure( Zinit[,sample.int(ncol(Zinit),m0*nChain,replace=TRUE)], dim=c(nParameter, m0, nChain), dimnames=list(parms=rownames(Zinit),NULL,NULL) )
            }
            # take the last row (or alternative samples with finite logDensity) and provide as x0 and logDensityComponents 0 
            resInitial <- .ensureFiniteLogDensityInInitialState(object, ZinitChains)
            sampleLogs <- do.call( .createInitialSampleLogs, c(list(object, sampleDimensions=sampleDimensions), resInitial) )
            ##value<< a SampleLogsObject with 1 sample and an initial population 
            sampleLogs
        })


#if(!exists(".computeInitialPopulation")) setGeneric(".computeInitialPopulation", function(object,...) standardGeneric(".computeInitialPopulation"))
#setMethod(".computeInitialPopulation", signature="PopulationSampler",
#' @export
computeInitialPopulation <-
        function( #object,
                ### Generate an initial population of states for \code{\link{twDEMCBlockInt}}.
                thetaPrior		##<< numeric vector (nParm) of point estimate of the mean
                ,covarTheta 	##<< numeric matrix (nParm x nParm) the covariance of parameters.
                ##<< Alternatively, can be given as vector (nParm) of the diagonal, i.e. variances, if all dimensions are independent
                ,nChainInPopulation = 1L   ##<< number of chains in population (defaults to nChain, i.e. one population)
                ,nPopulation = 1L       ##<< defaults to 1
                ,parameterNames = names(thetaPrior)
                ,m0= computeNRepresentativeRows(length(thetaPrior),nChainInPopulation)/(m0FiniteFac)
                ### number of initial states for each chain
                ,m0FiniteFac=1	##<< use a factor smaller than 1 to increase default m0 to account for only a portion of proposal results in finite densities
                ,doIncludePrior=TRUE
                ### If TRUE, then set last sample of chain 1 to the prior estimate, 
                ### which might be already a kind of best estimates by an optimization.
                ,subSpace=new("SubSpace")                       ##<< constraint of sampling range
        ){
                    ##<< scalar integer: number of chains, i.e. columns to generate
            nChain <- nPopulation * nChainInPopulation
            nParameter <- length(thetaPrior)
            Zinit <- if( nParameter==1L ){
                        abind( lapply( 1:nChain, function(i){ 
                                            samples <- rnorm( m0, mean=thetaPrior, sd=covarTheta)
                                            names(samples) <- parameterNames
                                            samples <- .letSubSpaceReflectBoundaries(subSpace, samples)
                                            matrix(samples, nrow=1, dimnames=list(parameterNames,NULL) ) })
                                , along=3 )
                    } else {
                        if( is.matrix(covarTheta))
                            abind( lapply( 1:nChain, function(i){ 
                                                samples <- rmvnorm( m0, mean=thetaPrior, sigma=covarTheta ) 
                                                colnames(samples) <- names(thetaPrior)
                                                samplesConstr <- .letSubSpaceReflectBoundariesMatrix(subSpace, samples)
                                                t(samplesConstr) 
                                            })
                                    , along=3 )
                        else
                            # independent dimensions, can sample each from rnorm
                            abind( lapply( 1:nChain, function(i){
                                                #iComp=1
                                                samples <- sapply( seq_along(thetaPrior), function(iComp){
                                                            rnorm(m0, mean=thetaPrior[iComp], sd=covarTheta[iComp])
                                                        })
                                                colnames(samples) <- names(thetaPrior)
                                                samplesConstr <- .letSubSpaceReflectBoundariesMatrix(subSpace, samples)
                                                t(samplesConstr)
                                            }), along=3 )
                    }
            # Set the last sample of chain 1 to the prior estimate, which might be already a kind of best estimates by an optimization.
            if( doIncludePrior )
                Zinit[,m0,1 ] <- .letSubSpaceReflectBoundaries(subSpace, thetaPrior)
            dimnames(Zinit) = list(parameterNames,NULL,NULL)
            Zinit
        }

if(!exists(".ensureFiniteLogDensityInInitialState")) setGeneric(".ensureFiniteLogDensityInInitialState", function(object,...) standardGeneric(".ensureFiniteLogDensityInInitialState"))
setMethod(".ensureFiniteLogDensityInInitialState", signature="PopulationSampler",
        function(object, Zinit){ 
            # compute logDensities of last state
            x0 <- adrop(Zinit[,ncol(Zinit), ,drop=FALSE],2)
            x0Chain <- x0[,1L]
            chainState <- computeChainState(object@chainSampler, x0Chain)   # invoke computation once outside alply for debugging
            logDensityComponents0 <- abind( alply(x0, 2, function(x0Chain){
                        chainState <- computeChainState(object@chainSampler, x0Chain)
                        getLogDensityComponents(chainState)
                    }), along=2)
            if( !any(is.finite(colSums(logDensityComponents0)))) stopDemac("All initial states yield non-finite logDensities. Check your logDensity functions.")
            if( getNLogDensityComponent(getBlockDimensions(getChainSampler(object))) != 0L ){
                #logDensityComponents0[1,1:3] <- -Inf
                iNonFinite <- which( !apply(logDensityComponents0,2, function(logDensityComponents0Chain){
                                    all(is.finite(logDensityComponents0Chain))}))
                if( length(iNonFinite) ){
                    # take n alternative samples from other rows and calculte logDensity
                    # then replace the non-Finite cases in x0
                    propNonfinite <- length(iNonFinite)/ncol(logDensityComponents0)
                    nAlternativeSamples <- ceiling(length(iNonFinite)*3/(1-propNonfinite))  
                    iSteps <- sample.int(ncol(Zinit)-1, nAlternativeSamples, replace=TRUE)
                    iChains <- sample.int( dim(Zinit)[3], nAlternativeSamples, replace=TRUE)
                    x0Alt <- abind( lapply(1:nAlternativeSamples, function(i){Zinit[,iSteps[i],iChains[i]]}), along=2)
                    logDensityComponents0Alt <- abind( alply(x0Alt, 2, function(x0Chain){
                                        chainState <- computeChainState(object@chainSampler, x0Chain)
                                        getLogDensityComponents(chainState)
                            }),along=2)
                    iFiniteAlt <- which( apply(logDensityComponents0Alt,2, function(logDensityComponents0Chain){
                                        all(is.finite(logDensityComponents0Chain))}))
                    if( length(iFiniteAlt) < length(iNonFinite) ){
                        #recover()
                        stop("could not get enough initial states with finite logDensity.")                        
                    } 
                    iFiniteAlt <- iFiniteAlt[1:length(iNonFinite)]
                    x0[,iNonFinite] <- x0Alt[,iFiniteAlt]
                    logDensityComponents0[,iNonFinite] <- logDensityComponents0Alt[,iFiniteAlt]
                }
            }
            ##value<< list with entries
            list(
                    initialPopulation=Zinit[,-ncol(Zinit), ,drop=FALSE]  ##<< array of parameters (nParameter, m0, nChain)
                    ,x0=x0                                               ##<< vector of parameters (nParameter)
                    ,logDensityComponents0=logDensityComponents0         ##<< logDensityComponent of x0
            )
        })

if(!exists(".createInitialSampleLogs")) setGeneric(".createInitialSampleLogs", function(object,...) standardGeneric(".createInitialSampleLogs"))
setMethod(".createInitialSampleLogs", signature="PopulationSampler", 
        function( object
                ### Generate an initial population of states for \code{\link{twDEMCBlockInt}}.
                ,sampleDimensions       
                ,initialPopulation
                ,x0
                ,logDensityComponents0
        ){
            sampleLogs <- object@sampleLogs     # template object
            sampleDimensions(sampleLogs) <- sampleDimensions
            initialPopulation(sampleLogs) <- initialPopulation
            sampleLogs <- setSampleRecord(sampleLogs,
                    iSample = 1L
                    ,parameters = x0
                    ,logDensityComponents = logDensityComponents0
            )
            sampleLogs            
        })

if(!exists("sampleLogs<-")) setGeneric("sampleLogs<-", function(object,value) standardGeneric("sampleLogs<-"))
# need to export, because of SannSampler inherits .adjustThinningToSpectralDensity from BatchSampler
#' @export
setReplaceMethod("sampleLogs", signature=c("PopulationSampler", "SampleLogs"), 
        function(object, value) {
            if( !(.isSampleLogsConsistentWithChainSampler(value, object@chainSampler)) ) 
                stop("Can only set sampleLogs that is consistent with ChainSampler.")
            sampleLogsOrig <- object@sampleLogs 
            object@sampleLogs <- value
            blockDimensions(object@sampleLogs) <- getBlockDimensions(object@chainSampler)
            hasChangedDimensions <-
                    getNParameter(value) != getNParameter(sampleLogsOrig) || 
                    any(getParameterNames(value) != getParameterNames(sampleLogsOrig)) ||
                    getNPopulation(value) != getNPopulation(sampleLogsOrig)                                         
            if( hasChangedDimensions ){
                if( getNParameter(sampleLogsOrig) )
                    warning("sampleLogs<-: replacing subSpacePopulations by empty objects because of changed dimensions")
                # changes sampleDimensions, need to reinitialized indexBounds
                # NULL objects for subspaces that should be set separately
                # setter takes care of initializing indexedBounds
                subSpacePopulations(object) <- rep(list(new("SubSpace")) 
                        , getNPopulation(getSampleDimensions(object)))
                object@nRepresentativeRows <- computeNRepresentativeRows(getNParameter(object@sampleLogs), getNChainInPopulation(object@sampleLogs))
                # temperatures and jumpProposer will be initialized in setupAndSample
            }
            object        
        })

.isSampleLogsConsistentWithChainSampler <- function(sampleLogs, chainSampler){
    isBlockDimensionsConsistentWith(getBlockDimensions(chainSampler), sampleLogs )
}

if(!exists("sampleLimiting")) setGeneric("sampleLimiting", function(object,...) standardGeneric("sampleLimiting"))
#' @export
setMethod("sampleLimiting", signature=c("PopulationSampler"),
        function(
                ### continue sampling aber burnin
                object                                  ##<< the Population sampler object
                ,nSample=3L*nRepresentativeRows         ##<< integer vector (nPopulation): number of samples to generate
                ,nSampleInitial=2*nRepresentativeRows   ##<< integer vector (nPopulation): number of samples to keep from end of previous sample
                ,nRepresentativeRows = computeNRepresentativeRows(getNParameterWithProposal(object@sampleLogs), getNChainInPopulation(object@sampleLogs))    
        ){
            sampleLogs(object) <- subsetTail( getSampleLogs(object), nSample=nSampleInitial)
            isBurnin(object) <- FALSE
            object <- setupAndSample(object, nSample=nSample)
            object
        })

#trace(setupAndSample, recover, signature = c("PopulationSampler"))    #untrace("setupAndSample", signature = c("PopulationSampler"), where =getNamespace("blockDemac"))
if(!exists("setupAndSample")) setGeneric("setupAndSample", function(object,...) standardGeneric("setupAndSample"))
#' @export
setMethod("setupAndSample", signature=c("PopulationSampler"),
        function(
                ### setup and perform the sampling
                object          ##<< the Population sampler object
                ,nSample        ##<< integer vector (nPopulation): number of samples to generate
                ,thin=4L        ##<< integer scalar: thinning interval used, i.e. after how many generations, a sample is recorded
                ,TStart = getEndTemperaturesPopulations(object@sampleLogs)     ##<< numeric scalar: initial temperature (or vector for components, or matrix for different populations)
                ,TEnd = TStart  ##<< numeric scalar: end temperature (or vector for components, or matrix for different populations)
                ,nIntervalPerRange = as.integer(round(computeNRepresentativeRows(getNParameter(object@sampleLogs), getNChainInPopulation(object@sampleLogs))/8))    ##<< combine sampling jumps for the as many thinning intervals
        ){
            force(TStart)   # else Start temperature is calcualted for extended SamplLogs
            if( length(nSample)==1 ) nSample <- rep(nSample, getNPopulation(object@sampleLogs))
            # initialize thinning interval of chain sampler 
            thin(object@chainSampler) <- thin
            #subSpacePopulations(object) <- subSpacePopulations
            object@nSampleBeforeSampling <- getNSamplePopulations(object@sampleLogs)
            object@nSampleBatchPopulations <- nSample
            object@nIntervalPerRange <- nIntervalPerRange
            # allocate space for new samples
            nSamplePopulations(object@sampleLogs) <- object@nSampleBeforeSampling + object@nSampleBatchPopulations
            # setup step temperatures
            object <- .initializeStepTemperatures(object, TStart=TStart,TEnd=TEnd)
            # set chain states to end of sampleLog
            x0 <- getParametersForPopulation(object@sampleLogs,1L)[,object@nSampleBeforeSampling[1L],1L]
            object@chainStates <- .initChainStates(
                  object@sampleLogs
                , iSamplePopulations=object@nSampleBeforeSampling
                , intermediateIds = getIntermediateIdsUsed(object@chainSampler)
                )
            object@acceptanceTracker <- resetAcceptanceTracker(object@acceptanceTracker
                    ,sampleDimensions=getSampleDimensions(object@sampleLogs)
                    , nSample=object@nSampleBatchPopulations
                    , nIntervalPerRange=object@nIntervalPerRange)
            object <- .resetJumpProposer(object, isBurnin=object@isBurnin)    # depends on thinning 
            # export things that do not change over sampling to remote processes
            if (.isParallel(object@cluster)){
                mcCommon <- .getMcCommon(object)
                #clusterExport(object@cluster, c("mcCommon", ".sampleRangeOneChain"), environment())  
                clusterExport(object@cluster, c("mcCommon"), environment())  
            } 
            object <- .samplePopulations(object)
            object
        })


#trace(blockDemac:::.initializeStepTemperatures, recover, signature = c("PopulationSampler"), where =getNamespace("blockDemac"))    #untrace("blockDemac:::.initializeStepTemperatures", signature = c("PopulationSampler"), where =getNamespace("blockDemac"))
if(!exists(".initializeStepTemperatures")) setGeneric(".initializeStepTemperatures", function(object,...) standardGeneric(".initializeStepTemperatures"))
setMethod(".initializeStepTemperatures", signature=c("PopulationSampler"),
        function( object
                , TStart = 1
                , TEnd = 1
        ){  
            object@stepTemperaturesPopulations <- .computeTemperatureDecreasePopulations(object, TStart, TEnd)
            thin <- getThin(object@chainSampler)
            thinnedStepTemperatures <- lapply(object@stepTemperaturesPopulations, function(stepTemperatures){
                        iSteps <- seq(thin,ncol(stepTemperatures),by=thin)
                        stepTemperatures[,iSteps ,drop=FALSE]
                    })
            stepTemperaturesPopulations(object@sampleLogs, iSample0=object@nSampleBeforeSampling) <- thinnedStepTemperatures             
            object
        })

if(!exists(".computeTemperatureDecreasePopulations")) setGeneric(".computeTemperatureDecreasePopulations", function(object,...) standardGeneric(".computeTemperatureDecreasePopulations"))
setMethod(".computeTemperatureDecreasePopulations", signature=c("PopulationSampler"),
        function( object
                ### Calculates the temperature for an exponential decrease from \code{TStart} to \code{Tend} after \code{nGen} steps. 	
                , TStart=1	    ##<< numeric matrix (nLogDensityComponents, nPopulation): the initial temperature (before the first step at iGen=0)
                , TEnd=1	##<< numeric matrix (nLogDensityComponents, nPopulation): the initial temperature (before the first step at iGen=0)
        ){
            if( !length(object@nSampleBatchPopulations) ) stop("nSampleBatchPopulations not initialized.")
            nGen <- object@nSampleBatchPopulations * getThin(object@chainSampler)
            nPop <- length(nGen)
            nLogDensityComponent <- getNLogDensityComponent(object@sampleLogs)
            if( nLogDensityComponent == 0){
                # return matrices with 0 rows
                res <- lapply(1:nPop, function(iPop){
                            matrix( NA_real_, nrow=0, ncol=nGen[iPop])
                        })
                return( res )
            }
            TStartM <- .formatPopulationTemperatures(object, TStart)
            TEndM <- .formatPopulationTemperatures(object, TEnd)
            lapply( 1:nPop, function(iPop){
                        t(.computeExponentialDecreaseVector(TStart=TStartM[,iPop],TEnd=TEndM[,iPop],nGen=nGen[iPop]))  
                    })
        })

if(!exists(".formatPopulationTemperatures")) setGeneric(".formatPopulationTemperatures", function(object,...) standardGeneric(".formatPopulationTemperatures"))
setMethod(".formatPopulationTemperatures", signature=c("PopulationSampler"),
        function(object,
                ### convert temperature specification to proper format for each population and logDensity component
                temperature     ##<< temperature (matrix, vector, scalar)
        ){
            ans <- temperature
            nPop <- getNPopulation(getSampleLogs(object))
            logDensityComponentNames <- getLogDensityComponentNames(getSampleLogs(object))
            nLogDensityComponent <- length(logDensityComponentNames)
            ##details<< if TStart or TEnd are given as vector, assume they refer to all populations.
            if( !is.matrix(ans) ) ans <- matrix(ans, nrow=length(ans), ncol=nPop)
            if( nLogDensityComponent > 1){
                ##details<< if rows in TStart or TEnd is one, assume they refer to all logDensityComponents
                if( nrow(ans)==1 ) ans <- matrix( ans, byrow=TRUE, nrow=nLogDensityComponent, ncol=ncol(ans))
            }
            if( nrow(ans) != nLogDensityComponent ) stop("number of rows in temperature must correspond to number of LogDensityComponents.")
            if( ncol(ans) != nPop ) stop("number of columns in temperature must correspond to number of populations.")
            rownames(ans) <- logDensityComponentNames
            ans
        })


if(!exists(".resetJumpProposer")) setGeneric(".resetJumpProposer", function(object,...) standardGeneric(".resetJumpProposer"))
setMethod(".resetJumpProposer", signature="PopulationSampler", 
        function(object, isBurnin=TRUE) {
            iParametersThatRequireJumpProposal <- getIParametersThatRequireJumpProposal(object@chainSampler)
            object@jumpProposer <- initializeDimensionalSettings(object@jumpProposer
                    , sampleDimensions=getSampleDimensions(object)
                    , thin=getThin(object@chainSampler)
                    , iParametersThatRequireJumpProposal=iParametersThatRequireJumpProposal )
            object@jumpProposer@isBurnin <- isBurnin
            object
        })



if(!exists(".samplePopulations")) setGeneric(".samplePopulations", function(object, ...) standardGeneric(".samplePopulations"))
setMethod(".samplePopulations", signature=c("PopulationSampler"),
        function(object){
            # (each nIntervalPersRange thinning intervals times thin generations)
            # Instead of handling all populations separately, in each step all participating  
            # chains are combined to a big array, to best use parallel calculation.
            # However, some populations may need fewer thinning steps. 
            # Hence, track is kept of the chains that still take part.
            iSampleInBatch <- 1L
            # populations that take part in this step (some may need fewer steps)
            iPopulationsStep <- which( object@nSampleBatchPopulations >= iSampleInBatch)   
            while( length(iPopulationsStep) ){    # as long as some populations participating (all thinning steps)
                # thinning interval after which fewer populations take part
                nSampleBatchMin <- min( object@nSampleBatchPopulations[iPopulationsStep] )   
                while( iSampleInBatch <= nSampleBatchMin ){
                    # anticipate that nIntervalPerRange may overshooting nSampleBatchMin
                    nThinningIntervalStep <- min(object@nIntervalPerRange, nSampleBatchMin - iSampleInBatch + 1L )
                    .signalProgress(object@progressOutput)
                    object <- .proposeSampleAndLogRange(object
                            ,iPopulationsStep=iPopulationsStep
                            ,iSampleInBatch=iSampleInBatch
                            ,nThinningIntervalStep=nThinningIntervalStep
                    )
                    iSampleInBatch <- iSampleInBatch + nThinningIntervalStep        
                } # end not change in participating populations
                # update populations that take part in the following sequence of steps
                iPopulationsStep <- which( object@nSampleBatchPopulations >= iSampleInBatch)
                #if( length(iPopulationsStep) == 1 ) recover()        
            } # end main loop
            .signalProgress(object@progressOutput, appendLF = TRUE)
            object
        }
)

.tmp.f <- function(){
    #library(profr)
    p1 <- profr({
                for( i in 1:10 ){
                .proposeSampleAndLogRange(object
                        ,iPopulationsStep=iPopulationsStep
                        ,iSampleInBatch=iSampleInBatch
                        ,nThinningIntervalStep=nThinningIntervalStep
                )            
                }
            })
    plot(p1)
}

.signalProgress <- function(progressOutput, appendLF = FALSE){
    if( length(progressOutput) ) message(progressOutput, appendLF=appendLF)  
}

if(!exists(".getMcCommon")) setGeneric(".getMcCommon", function(object, ...) standardGeneric(".getMcCommon"))
setMethod(".getMcCommon", signature=c("PopulationSampler"),
        function(object){
            mcCommon <- list(
                    chainSampler = object@chainSampler
                    ,subSpacePopulations = object@subSpacePopulations
                    ,temperatures = object@stepTemperaturesPopulations
                    ,sampleDimensions = getSampleDimensions(object@sampleLogs)
            )
        })
        

if(!exists(".proposeSampleAndLogRange")) setGeneric(".proposeSampleAndLogRange", function(object,...) standardGeneric(".proposeSampleAndLogRange"))
setMethod(".proposeSampleAndLogRange", signature=c("PopulationSampler"),
        function(object, iPopulationsStep, iSampleInBatch, nThinningIntervalStep){  # as long as the same populations take part
            iChainsStep <- getIChainsForPopulations(getSampleDimensions(object@sampleLogs), iPopulationsStep)
            pAcceptPops  <- computePopulationAcceptanceRates(object@acceptanceTracker)
            iPopulationChains <- getIPopulationsForChains(getSampleDimensions(object@sampleLogs), iChainsStep)
            pAcceptChains <- pAcceptPops[iPopulationChains]
            # jumps relate already only to the populations taking part in this range (iPopulationsStep)
            object@jumpProposer <- adjustToAcceptanceRate(object@jumpProposer, acceptanceTracker = object@acceptanceTracker)
            jumps <- proposeJumps( object@jumpProposer
                    , nGeneration=nThinningIntervalStep*getThin(object@chainSampler)
                    , sampleLogs=object@sampleLogs
                    , iPopulationsStep = iPopulationsStep
                    , iCurrentSamples = object@nSampleBeforeSampling+iSampleInBatch-1  # current state before the current index
            )
            # samplesChains also takes information of the participating populations
            samplesChains <- sampleRangeAllChains(object
            #samplesChains <- blockDemac:::sampleRangeAllChains(object
                    , chainStates = object@chainStates[iChainsStep]
                    , intervalInfoChains =list(
                            iChainsStep=iChainsStep
                            ,iSample0=iSampleInBatch-1
                            ,nSample=nThinningIntervalStep
                            ,step=jumps$jump
                            , rExtra=jumps$logSnookerDensityMultiplier
                            , pAccept=pAcceptChains
                    ), mcCommon=blockDemac:::.getMcCommon(object))
            # samplesChains holds only the participating chains
            object@chainStates[iChainsStep] <- lapply( samplesChains, "[[", "chainState" )
            sampleLogChains <- lapply( samplesChains, "[[", "sampleLog" )
            object@sampleLogs <- recordSampleLogChains(object@sampleLogs, sampleLogChains
                    , iPopulations=iPopulationsStep, iSample0=object@nSampleBeforeSampling+iSampleInBatch-1 )
            proportionAcceptedInInterval <- abind( lapply(sampleLogChains, getProportionAcceptedInInterval), rev.along=0)
            object@acceptanceTracker <- recordAcceptanceRate(object@acceptanceTracker, proportionAcceptedInInterval=proportionAcceptedInInterval, iChains=iChainsStep)
            object
        })


.initChainStates <- function(
        ### Initialize chainState to given recored sample
        sampleLogs         ##<< SampleLogs Object
        ,iPopulations=1:getNPopulation(sampleLogs)      ##<< integer vector of populations that take part 
        ,iSamplePopulations=getNSamplePopulations(sampleLogs)[iPopulations]  ##<< integer vector (length iPopulations): index of the sample to put to chainState for each population
        ,intermediateIds=character(0)   ##<< names of the intermediates, that will be initialized to "not up to date, i.e. empty list"         
        ,chainStates         ##<< list of ChainState objects to initialize
){
    if( missing(chainStates) ){
        chainStates <- rep(list(new("ChainState",getParametersAndInitialForPopulation(sampleLogs,1L)[,1L,1L] )), getNChain(sampleLogs))
    }
    iChains <- as.vector(sapply( iPopulations, function(iPop){ getIChainsForPopulation(sampleLogs,iPop) }))
    sampleLogChains <- getSampleLogChainsForPopulations(sampleLogs, iPopulations)
    iSampleChains <- rep(iSamplePopulations,each=getNChainInPopulation(sampleLogs))
    intermediates <- structure( rep( list(list()), times=length(intermediateIds)), names=intermediateIds )
    updatedChainStates <- mapply( function(chainState, sampleLogChain, iSample){
                intermediates(chainState) <- intermediates
                initializeFromSampleLogChain(chainState, sampleLogChain, iSample)
            }, chainState=chainStates[iChains], sampleLogChain=sampleLogChains, iSample=iSampleChains)
    updatedChainStates
    
}

if(!exists("sampleBurnedIn")) 


if(!exists("getSampleDimensions")) setGeneric("getSampleDimensions", function(object) standardGeneric("getSampleDimensions"))
setMethod("getSampleDimensions", signature="PopulationSampler", function(object) {getSampleDimensions(object@sampleLogs)})


modifyIntermediateInSampler <- function(
    ### change something in intermediate of given population sampler
    sampler             ##<< PopulationSampler object
    ,intermediateName   ##<< name of the intermediate to change
    , fModify=function(intermediate,...){intermediate}  ##<< function that receives an IntermediateUpdater object and returns one
    ,...    ##<< further arguments to fModify
){
    sampler@chainSampler@blockUpdaters@intermediateUpdaters@intermediateUpdaters[[intermediateName]] <-
            fModify(sampler@chainSampler@blockUpdaters@intermediateUpdaters@intermediateUpdaters[[intermediateName]])
    ##value<< a modified sampler object keeping all sampleLogs
    sampler
}
attr(modifyIntermediateInSampler,"modifyIntermediateInSampler") <- function(){
    if( FALSE ) {
        # sometimes an update of the additional parameters to the intermediate is required before continued sampling
        sampler2 <- modifyIntermediateInSampler(sampler, "deltaA_rich", function(im){
                    im@argsUpdateFunction$minSd2DiscrFractionVal <- 1/20
                    im
                })
    }
}

if(!exists("setNoCluster")) setGeneric("setNoCluster", function(object,...) standardGeneric("setNoCluster"))
#' @export
setMethod("setNoCluster", signature="PopulationSampler", 
        function( object
        ### Set to not usiung any cluster
        ){
            object@cluster <- structure(list(),class=c("SOCKcluster","cluster"))
            object
        })


if(!exists("getBlockDimensions")) setGeneric("getBlockDimensions", function(object) standardGeneric("getBlockDimensions"))
setMethod("getBlockDimensions", signature="PopulationSampler", function(object) { getBlockDimensions(object@chainSampler)})


#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("PopulationSampler")
if(!exists("getChainStates")) setGeneric("getChainStates", function(object) standardGeneric("getChainStates"))
#' @export
setMethod("getChainStates", signature="PopulationSampler", function(object) {object@chainStates})
#if(!exists("chainStates<-")) setGeneric("chainStates<-", function(object,value) standardGeneric("chainStates<-"))
#setReplaceMethod("chainStates", signature=c("PopulationSampler", "list"), function(object, value) {object@chainStates <- value; object})

if(!exists("getSampleLogs")) setGeneric("getSampleLogs", function(object) standardGeneric("getSampleLogs"))
#' @export
setMethod("getSampleLogs", signature="PopulationSampler", function(object) {object@sampleLogs})

if(!exists("getChainSampler")) setGeneric("getChainSampler", function(object) standardGeneric("getChainSampler"))
#' @export
setMethod("getChainSampler", signature="PopulationSampler", function(object) {object@chainSampler})
if(!exists("chainSampler<-")) setGeneric("chainSampler<-", function(object,value) standardGeneric("chainSampler<-"))
#' @export
setReplaceMethod("chainSampler", signature=c("PopulationSampler", "ChainSamplerImpl"), function(object, value) {
            if( !(.isSampleLogsConsistentWithChainSampler(object@sampleLogs, value)) ){ 
                warning("Setting inconsistent ChainSampler invalidates sampleLogs. Need to set sampleLogs again.")
                object@sampleLogs <- new("SampleLogs")
            }
            object@chainSampler <- value
            blockDimensions(object@sampleLogs) <- getBlockDimensions(object@chainSampler)
            object            
        })

if(!exists("jumpProposer<-")) setGeneric("jumpProposer<-", function(object,value) standardGeneric("jumpProposer<-"))
#' @export
setReplaceMethod("jumpProposer", signature=c("PopulationSampler", "JumpProposer"), function(object, value) {
            object@jumpProposer <- value
            object <- .resetJumpProposer(object)
            object
        })

if(!exists("getAcceptanceTracker")) setGeneric("getAcceptanceTracker", function(object) standardGeneric("getAcceptanceTracker"))
#' @export
setMethod("getAcceptanceTracker", signature="PopulationSampler", function(object) {object@acceptanceTracker})
if(!exists("acceptanceTracker<-")) setGeneric("acceptanceTracker<-", function(object,value) standardGeneric("acceptanceTracker<-"))
#' @export
setReplaceMethod("acceptanceTracker", signature=c("PopulationSampler", "AcceptanceTracker"), function(object, value) {
            stop("not implemented yet. initialize object")
            object@acceptanceTracker <- value; object
        })

if(!exists("getProgressOutput")) setGeneric("getProgressOutput", function(object) standardGeneric("getProgressOutput"))
#' @export
setMethod("getProgressOutput", signature="PopulationSampler", function(object) {object@progressOutput})
if(!exists("progressOutput<-")) setGeneric("progressOutput<-", function(object,value) standardGeneric("progressOutput<-"))
#' @export
setReplaceMethod("progressOutput", signature=c("PopulationSampler", "character"), function(object, value) {object@progressOutput <- value; object})

if(!exists("getSubSpacePopulations")) setGeneric("getSubSpacePopulations", function(object) standardGeneric("getSubSpacePopulations"))
setMethod("getSubSpacePopulations", signature="PopulationSampler", function(object) {object@subSpacePopulations})
if(!exists("subSpacePopulations<-")) setGeneric("subSpacePopulations<-", function(object,value) standardGeneric("subSpacePopulations<-"))
setReplaceMethod("subSpacePopulations", signature=c("PopulationSampler", "list"), function(object, value) {
            if( length(value) != getNPopulation(object@sampleLogs)) stop("Need to provide a subspace for each population.")
            if( !all(sapply(value, is, "SubSpace"))) stop("Need to provide SubSpace objects.")
            subSpaces <- lapply( value, initializeIndexedBoundsForParameterNames, parameterNames=getParameterNames(object@sampleLogs))            
            object@subSpacePopulations <- subSpaces 
            object
        })

if(!exists("isBurnin")) setGeneric("isBurnin", function(object,...) standardGeneric("isBurnin"))
#' @export
setMethod("isBurnin", signature="PopulationSampler", function(object,...
        ### Getter method for slot isBurnin
        ) {object@isBurnin})
if(!exists("isBurnin<-")) setGeneric("isBurnin<-", function(object,value) standardGeneric("isBurnin<-"))
#' @export
setReplaceMethod("isBurnin", signature=c("PopulationSampler", "logical"), function(object, value
        ### Setter method for slot isBurnin
        ) {object@isBurnin <- value; object})

if(!exists("getNRepresentativeRows")) setGeneric("getNRepresentativeRows", function(object,...) standardGeneric("getNRepresentativeRows"))
#' @export
setMethod("getNRepresentativeRows", signature="PopulationSampler", function(object,...
        ### Getter method for slot nRepresentativeRows
        ) {object@nRepresentativeRows})




