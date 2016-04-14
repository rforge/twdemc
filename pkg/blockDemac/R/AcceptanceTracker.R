#' @include AcceptanceTrackerSettings.R
#' @include SampleDimensionsImpl.R
#' @include SampleLogChainImpl.R

setClass("AcceptanceTracker"
        ,representation(  
                settings="AcceptanceTrackerSettings"
                ,proportionAcceptedWindow="array"
                ,iAcceptWindow="integer"	##<< index within acceptance window
                ,acceptanceLog="array"	# 0..1: (nBlock, nSample, nChain)
                ,countSample="integer"	# number of recorded steps of acceptance rates
                ,sampleDimensions = "SampleDimensions"
        )
        ,prototype(sampleDimensions=new("SampleDimensionsImpl"))
)

if(!exists("resetAcceptanceTracker")) setGeneric("resetAcceptanceTracker", function(object,...) standardGeneric("resetAcceptanceTracker"))
#' @export
setMethod("resetAcceptanceTracker", signature="AcceptanceTracker", 
        function( object
                ### Reset Tracking to given sample dimensions. 
                ,sampleDimensions       ##<< class SampleDimensions: the dimensions of the sampleLog
                ,nSample                ##<< integer vector (nPopulation): number of samples in batch
                ,nIntervalPerRange      ##<< integer scalar: number of thinning interval process in one record
        ){
            object@sampleDimensions <- sampleDimensions
            object <- .initAcceptanceWindow(object, sampleDimensions, nIntervalPerRange)
            object <- .initializeAcceptanceLog(object, sampleDimensions, nSample)
            object
        }
)


if(!exists(".initAcceptanceWindow")) setGeneric(".initAcceptanceWindow", function(object,...) standardGeneric(".initAcceptanceWindow"))
setMethod(".initAcceptanceWindow", signature="AcceptanceTracker", 
        function( object
                ,sampleDimensions
                ,nIntervalPerRange
        ){
            nChain <- getNChain(sampleDimensions)
            # proportionAcceptedWindow holds the number of proportion of accepted steps for each thinning interval for each chain
            # in each thinning interval a value is added at the end
            # If the windows' boundaries are filled up
            # , the last nGenBack states are copied to the first interval
            # hence we need space for 2*acceptanceWindowWidth
            acceptWindowWidth <- as.integer(object@settings["widthAcceptWindow"]) 
            ##value<< numeric array (width x block x chain ), with number of accepted steps per thinning interval, 
            ## initialized with settings["initialAcceptanceRate"]
            nBlock <- getNBlock(sampleDimensions)
            blockNames <- getBlockNames(sampleDimensions)
            object@proportionAcceptedWindow <- array( NA_real_ 
                    , dim=c(nBlock, 2L*acceptWindowWidth+nIntervalPerRange,nChain ) 
                    , dimnames=list(blocks=blockNames, steps=NULL, chains=NULL) )
            object@iAcceptWindow <- 1L
            object@proportionAcceptedWindow[,object@iAcceptWindow,] <-  object@settings["initialAcceptanceRate"]
            object
        }
)

if(!exists(".initializeAcceptanceLog")) setGeneric(".initializeAcceptanceLog", function(object,...) standardGeneric(".initializeAcceptanceLog"))
setMethod(".initializeAcceptanceLog", signature="AcceptanceTracker", 
        function( object
                ,sampleDimensions
                ,nSample
        ){
            nChain <- getNChain(sampleDimensions)
            #nSample <- max(getNSamplePopulations(sampleDimensions))
            nBlock <- getNBlock(sampleDimensions)
            blockNames <- getBlockNames(sampleDimensions)
            object@acceptanceLog <- array( NA_real_, c(nBlock, max(nSample), nChain)
                    ,dimnames=list(blockNames, steps=NULL, chains=NULL))
            object@countSample <- 0L
            object@acceptanceLog[,object@countSample,] <- object@settings["initialAcceptanceRate"]
            object
        }
)


if(!exists("recordAcceptanceRate")) setGeneric("recordAcceptanceRate", function(object,...) standardGeneric("recordAcceptanceRate"))
setMethod("recordAcceptanceRate", signature=c("AcceptanceTracker"),
        function(object
                ### recored sumAccepted for steps iThin and successive length(sumAccepted) thinning intervals
                ,proportionAcceptedInInterval    ##<< array (nBlock, nInterval, nChainRange) of sum of acceptance over each thinning interval
                ,iChains        ##<< index of chains in sumAccepted
        ){
            if( missing(proportionAcceptedInInterval)) stop("Need to provide proportionAcceptedInInterval.")
            if( missing(iChains) ) iChains <- 1:getNChain(object@sampleDimensions)
            nInterval <- ncol(proportionAcceptedInInterval)
            assert_that( object@countSample + nInterval <= ncol(object@acceptanceLog) )
            acceptWindowWidth <- as.integer(object@settings["widthAcceptWindow"]) 
            object@proportionAcceptedWindow[,object@iAcceptWindow+(1:nInterval), iChains] <- proportionAcceptedInInterval
            for( i in 1:nInterval){
                step <- object@iAcceptWindow + i
                nThinBack <- min(acceptWindowWidth,step)     
                curWindow <- step+1-(nThinBack:1)
                object@acceptanceLog[,object@countSample+i,iChains] <- pAcceptChains <- 
                        apply(object@proportionAcceptedWindow[,curWindow,iChains ,drop=FALSE],c(1,3),sum) /
                        (nThinBack)
            }
            object@iAcceptWindow <- object@iAcceptWindow + nInterval
            object@countSample <- object@countSample + nInterval
            #        
            #row in acceptance Window to record acceptance, if reaches or exceeds window length, then 
            # copy second part to first (rewind)
            nExceed <- object@iAcceptWindow - (2L*acceptWindowWidth) #
            if( nExceed >= 0){	
                # works in R: part of destination is used as source
                .tmp.old <- function(){
                    object@proportionAcceptedWindow[,1:(acceptWindowWidth+nExceed),] <- object@proportionAcceptedWindow[ 
                            ,acceptWindowWidth+(1:(acceptWindowWidth+nExceed)),]
                    object@proportionAcceptedWindow[,(acceptWindowWidth+nExceed+1):ncol(object@proportionAcceptedWindow),] <- NA
                    object@iAcceptWindow <- (acceptWindowWidth+nExceed)
                }
                object@proportionAcceptedWindow[,1:acceptWindowWidth,] <- object@proportionAcceptedWindow[ 
                        ,acceptWindowWidth+nExceed +(1:acceptWindowWidth),]
                object@proportionAcceptedWindow[,(acceptWindowWidth+1):ncol(object@proportionAcceptedWindow),] <- NA
                object@iAcceptWindow <- acceptWindowWidth                
            }
            object
        }
)

# note the "missing" iPopulations in the signature here compared with the method below
if(!exists("computePopulationAcceptanceRates")) setGeneric("computePopulationAcceptanceRates", function(object,iPopulations,...) standardGeneric("computePopulationAcceptanceRates"))
#' @export
setMethod("computePopulationAcceptanceRates", signature=c("AcceptanceTracker","missing"),
        function(object
                ### Average Acceptance rates within across chains.
                , iPopulations  ##<< missing by signature
        ){
            if( object@countSample == 0L) return( rep(object@settings["initialAcceptanceRate"],getNPopulation(object@sampleDimensions)))
            pAcceptChains <- adrop(object@acceptanceLog[,object@countSample, ,drop=FALSE],2) # last recorded row
            # acceptance rate of populations as a combination of average across chains  
            # and minimum of the blocks
            as.vector(popApplyTwDEMC(pAcceptChains, getNPopulation(object@sampleDimensions), .computeAcceptanceRateAcrossChainsAndGroups ))
        }
)       
#' @export
setMethod("computePopulationAcceptanceRates", signature=c("AcceptanceTracker","integer"),
        function(object
                ### Average Acceptance rates across chains within a population.
                , iPopulations  ##<< index of populations
        ){
            if( object@countSample == 0L) return( rep(object@settings["initialAcceptanceRate"],length(iPopulations)))
            iChains <- getIChainsForPopulations(object@sampleDimensions, iPopulations) 
            pAcceptChains <- adrop(object@acceptanceLog[,object@countSample,iChains ,drop=FALSE],2) # last recorded row
            as.vector(popApplyTwDEMC(pAcceptChains, length(iPopulations), .computeAcceptanceRateAcrossChainsAndGroups ))
        }
)

if(!exists("computeChainGroupAcceptanceRate")) setGeneric("computeChainGroupAcceptanceRate", function(object,iPopulation,iChains,...) standardGeneric("computeChainGroupAcceptanceRate"))
#' @export
setMethod("computeChainGroupAcceptanceRate", signature=c("AcceptanceTracker","integer","integer"),
        function(object
                ### Average Acceptance rates across specified chains within a population.
                , iPopulation
                , iChains
        ){
            if( object@countSample == 0L) return( rep(object@settings["initialAcceptanceRate"],length(iChains)))
            iChainsGlobal <- getIChainsForPopulation(object@sampleDimensions, iPopulation)[iChains] 
            pAcceptChains <- adrop(object@acceptanceLog[,object@countSample,iChainsGlobal ,drop=FALSE],2) # last recorded row
            .computeAcceptanceRateAcrossChainsAndGroups(pAcceptChains)
        }
)

.computeAcceptanceRateAcrossChainsAndGroups <- function(
        ### compute a scalar acceptance rate across blocks and chains
        pAcceptChains       ##<< numeric matrix (nBlock x nChain)
        ,pQuantBlocks=0.2   ##<< quantile to choose acceptance rate across blocks of one chain    
){
    ##details<<
    ## Sometimes, onle a few special blocks have low acceptance rate, while other parameters are fine.
    ## Hence, instead of taking the minimum across blocks, a low quantile is used.
    # acceptance rate of populations as a combination of average and minimum across chains  
    # and before low quantile across the blocks
    blockMin <- apply(pAcceptChains,2,quantile,pQuantBlocks, na.rm=TRUE) 
    (mean(blockMin) + min(blockMin))/2
}




#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("AcceptanceTracker")
if(!exists("getAcceptanceLog")) setGeneric("getAcceptanceLog", function(object,...) standardGeneric("getAcceptanceLog"))
#' @export
setMethod("getAcceptanceLog", signature="AcceptanceTracker", function(object,...
        ### Getter method for slot acceptanceLog
        ) {object@acceptanceLog})

if(!exists("getCountSample")) setGeneric("getCountSample", function(object,...) standardGeneric("getCountSample"))
#' @export
setMethod("getCountSample", signature="AcceptanceTracker", function(object,...
        ### Getter method for slot countSample
        ) {object@countSample})


