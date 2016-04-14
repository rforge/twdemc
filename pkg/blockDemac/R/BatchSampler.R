#' @include PopulationSampler.R

#' @export
newBatchSampler <- function(
        ### create a new BatchSampler based on given blockSpecifications
        ...             ##<< arguments to \code{\link{newPopulationSampler}}
        ,thin           ##<< initial thinning interval
        , convergenceControls = compileBatchConvergenceControls()
        , className="BatchSampler"
){
    sampler <- newPopulationSampler(..., className=className)
    sampler@convergenceControls = convergenceControls
    sampler@thin <- thin
    sampler            
}

#' @export
compileBatchConvergenceControls <- function(
        ### compile list of convergence control parameters
        #,maxRelTChange=0.05 	##<< if Temperature of the components changes less than specified value, report result
        #,minBaseTemperature=1e-3    ##<< if Temperature-1  gets below this temperature, report result
        maxGelmanDiag=1.3			##<< do not change Temperature, if variance between chains is too high, i.e. Gelman Diag is above this value
        , upperSpecVarRatio=20		##<< if proprotion of spectral Density to Variation is higher than this value, signal problems 
        , lowerSpecVarRatio=8       ##<< if proprotion of spectral Density to Variation is lower than this value, decrease thinning interval in next batch  
        , qSpecVarRatio=0.9         ##<< scalar numeric: quantile of distribution of spectral Density to Variance ratios that should be compared to upperSpecVarRatio 
        , maxThin=64                ##<< scalar positive integer: maximum thinning interval. If not 0 then thinning is increased in case of too high spectral density
        , thinIncreaseFactor=1.4    ##<< scalar numeric: factor to increase thinning interval when spectral density is too high 
        , qLogDiffDrift=0.8         ##<< scalar numeric: quantile in chiSquare distribution to determine maxDrift (differing logDensity)
){
    ##<< list with entries corresponding to the function arguments.
    list(
            #,maxRelTChange=maxRelTChange
            #,minBaseTemperature=minBaseTemperature
            maxGelmanDiag=maxGelmanDiag
            ,upperSpecVarRatio=upperSpecVarRatio
            ,lowerSpecVarRatio=lowerSpecVarRatio
            ,qSpecVarRatio=qSpecVarRatio
            ,maxThin=maxThin
            ,thinIncreaseFactor=thinIncreaseFactor
            ,qLogDiffDrift=qLogDiffDrift
    )
}


#' @export
setClass("BatchSampler",contains="PopulationSampler" 
        ,representation(
                thin = "integer"	            ##<< current thinning interval
                ,convergenceDiagnostics="list"  ##<< recently evaluated diagnostics
                ,batchLog="data.frame"	        ##<< logging information on batches 
                ,convergenceControls="list"	    ##<< list of parameters controlling the convergence
                ,cntConsequtiveMinorDecreaseBatches="integer"	##<< count of consequtive batches that fulfilled convergence criteria
        )
        ,prototype(
                thin = 2L
                , convergenceDiagnostics = list()
                , convergenceControls = compileBatchConvergenceControls()
                , cntConsequtiveMinorDecreaseBatches = 0L
        )
)
#str3(new("BatchSampler"))


if(!exists("sampleBatches")) setGeneric("sampleBatches", function(object,...) standardGeneric("sampleBatches"))
#' @export
setMethod("sampleBatches", signature="BatchSampler", function(object
                , nBatchMax=4L      ##<< maximum number of batches to sample      
                , nSampleBatch=     ##<< number of samples in batch
                        as.integer(round( {if(isBurnin(object)) 3 else 0.5}*getNRepresentativeRows(object)))
                , TStart=getEndTemperaturesPopulations(getSampleLogs(object))
                , TEnd=TStart
                , nSampleDiagnostics = as.integer(round(
                        max( 50/getNChainInPopulation(getSampleLogs(object))
                        ,2*getNRepresentativeRows(object))))  ##<< number of samples in tail to calculate diagnostics on
                , restartFilename = character(0)    ##<< if supplied a fileName, then object is saved before each batch as "sampler" in file <name>_<iBatch>.RData and <name>_end.RData
        ){
            # rows are not independent (especially when spectral density is high
            # therefore need to take longer part than getNRepresentativeRows(object) 
            ## to calcualte diagnostics on.
            object0 <- object
            if( nSampleDiagnostics < getNRepresentativeRows(object) ) warning(
                        "BatchSampler.sampleBatches: specified low nSampleDiagnostics. Should be at least nRepresentativeRows=",getNRepresentativeRows(object) )
            isConstantTemperature <- all( TStart==TEnd )
            if( !object@isBurnin && !isConstantTemperature ) warning("BatchSampler.sampleBatches: specified decreasing temperature but contrary also specified that its no burnin.")
            batchTemperatures <- .computeBatchTemperatures(object, TStart, TEnd, nBatchMax)
            message("sampleBatches: sampling ",nBatchMax," batches of ",nSampleBatch," samples")
            iBatch <- 1L
            object@cntConsequtiveMinorDecreaseBatches <- 0L
            TStartBatch <- .formatPopulationTemperatures(object, TStart)
            logsTail <- .getDiagnosticsLogsTail(object, nSampleDiagnostics )
            object@convergenceDiagnostics <- diagnostics <- 
                    computeConvergenceDiagnostics(object, logsTail )
            object <- .logAndMessageConvergenceRecord(object, iBatch=0L )
            while( iBatch <= nBatchMax ){
                # restart ?
                if( length(restartFilename)==1 ){
                    sampler <- object
                    save(sampler, file=paste(restartFilename,"_",iBatch,".RData",sep=""))
                }
                object <- .adjustThinningToSpectralDensity(object, logsTail)
                TEndPops <- adrop(batchTemperatures[,iBatch, ,drop=FALSE],2L)
                object <- setupAndSample(object, nSample=nSampleBatch, thin=object@thin, TStart=TStartBatch, TEnd=TEndPops)
                logsTail <- .getDiagnosticsLogsTail(object, nSampleDiagnostics )
                object@convergenceDiagnostics <- diagnostics <- 
                        computeConvergenceDiagnostics(object, logsTail )
                object <- .logAndMessageConvergenceRecord(object, iBatch=iBatch )
                #plot(asMcmc(logsTail),smooth=FALSE)
                #plot(asMcmc(object@sampleLogs),smooth=FALSE)
                object@cntConsequtiveMinorDecreaseBatches <- if( isConstantTemperature && .isConverged(object) ) 
                    object@cntConsequtiveMinorDecreaseBatches + 1L else 0L
                if( object@cntConsequtiveMinorDecreaseBatches == 2L){   
                    .showFinishEarlyMessage(object)
                    break
                }
                TStartBatch <- TEndPops
                iBatch <- iBatch + 1L
            }
            if( length(restartFilename)==1 ){
                sampler <- object
                save(sampler, file=paste(restartFilename,"_end.RData",sep=""))
            }
            object
        })
        
if(!exists(".computeBatchTemperatures")) setGeneric(".computeBatchTemperatures", function(object,...) standardGeneric(".computeBatchTemperatures"))
#' @export
setMethod(".computeBatchTemperatures", signature="BatchSampler", 
        function(object
                ### Getter method for slot convergenceDiagnostics
                , TStart        ##<< 
                , TEnd
                , nBatchMax        
        ) {
            nPop <- getNPopulation(getSampleLogs(object))
            nGen <- rep(nBatchMax, nPop)
            nLogDensityComponent <- getNLogDensityComponent(getSampleLogs(object)) 
            if( nLogDensityComponent == 0){
                # return matrices with 0 rows
                res <- lapply(1:nPop, function(iPop){
                            matrix( NA_real_, nrow=0, ncol=nGen[iPop])
                        })
                return( res )
            }
            TStartM <- .formatPopulationTemperatures(object, TStart)
            TEndM <- .formatPopulationTemperatures(object, TEnd)
            ##value<< numeric array (nLogDencomp x nBatch x nPop) of end termperatures at current batch
            abind( tmp <- lapply( 1:nPop, function(iPop){
                                t(.computeExponentialDecreaseVector(TStart=TStartM[,iPop],TEnd=TEndM[,iPop],nGen=nGen[iPop]))  
                            }), along=3L)
        })

if(!exists(".getDiagnosticsLogsTail")) setGeneric(".getDiagnosticsLogsTail", function(object,...) standardGeneric(".getDiagnosticsLogsTail"))
#' @export
setMethod(".getDiagnosticsLogsTail", signature="BatchSampler", function(object
                ### get the part of the sample logs that convergence diagnostics should be calculated on
                , nSampleDiagnosticsC=as.integer(round(1.5*getNRepresentativeRows(object)))  ##<< number of samples in tail to calculate diagnostics on
        ){
            logsTail <- if( object@isBurnin && min(getNSamplePopulations(object@sampleLogs)) > nSampleDiagnosticsC ){
                        subsetTail(object@sampleLogs, nSample = nSampleDiagnosticsC)
                    }else {
                        object@sampleLogs
                    }
        })

if(!exists(".logAndMessageConvergenceRecord")) setGeneric(".logAndMessageConvergenceRecord", function(object,...) standardGeneric(".logAndMessageConvergenceRecord"))
#' @export
setMethod(".logAndMessageConvergenceRecord", signature="BatchSampler", function(object
                ### append a row to batchLog and display message
                , iBatch
                , isDisplayMessage=TRUE
        ){
            diagnostics = object@convergenceDiagnostics            
            logRecord <- c( list( iBatch=iBatch, thin=object@thin
                            , maxGelmanDiagChains=max(diagnostics$gelmanDiagChains)
                            , upperSpecVarRatio=diagnostics$upperSpecVarRatio
                            , logDen=median(diagnostics$logDenBlocksStartEnd[2L,])
                    )
                    , diagnostics[c("isLogDenDrift","gelmanDiagPops")]
            )[c("iBatch","logDen","isLogDenDrift","gelmanDiagPops","maxGelmanDiagChains","upperSpecVarRatio","thin")]
            object@batchLog <- rbind(object@batchLog, as.data.frame(logRecord))
            rownames(object@batchLog) <- NULL
            if( isDisplayMessage )
                message(paste(capture.output(print(tail(object@batchLog,4))),collapse="\n"))
            object
        })

if(!exists(".adjustThinningToSpectralDensity")) setGeneric(".adjustThinningToSpectralDensity", function(object,...) standardGeneric(".adjustThinningToSpectralDensity"))
#' @export
setMethod(".adjustThinningToSpectralDensity", signature="BatchSampler", function(object
                ### compute and set thinning interval so that spectral density of thinned log is within bounds 
                , logsTail
                , isThinLogs=TRUE
        ){
            thinnedTail <- logsTail
            #plot(asMcmc(logsTail),smooth=FALSE)
            thinFac <- 1
            nRepresentativeRows <- getNRepresentativeRows(object)
            diagnostics <- getConvergenceDiagnostics(object)
            ctrl <- getConvergenceControls(object)
            # progressively thin logsTail and assess spectral density
            # if thinning interval is increased, squeeze also the entire sampleLog
            while( (min(getNSamplePopulations(getSampleLogs(object)))/thinFac > nRepresentativeRows) && 
                    (isHighAutoCorrelation <- is.finite(diagnostics$upperSpecVarRatio) && 
                        (diagnostics$upperSpecVarRatio > ctrl$upperSpecVarRatio)) 
                    ){
                thinIncreaseFactor <- min(ctrl$thinIncreaseFactor
                        , min(getNSamplePopulations(getSampleLogs(object)))*thinFac / nRepresentativeRows )
                thinFac <- thinFac*thinIncreaseFactor
                thinnedTail <- squeeze(thinnedTail,fraction=1/thinIncreaseFactor  )
                #getNSamplePopulations(thinnedTail)
                #plot(asMcmc(thinnedTail),smooth=FALSE)
                diagnostics <- computeConvergenceDiagnostics(object, thinnedTail )
            }
            # if spectral density is very low, decrease thinning interval
            if(
                    (object@thin > 2) &&
                    (isLowAutoCorrelation <- is.finite(diagnostics$upperSpecVarRatio) && 
                        (diagnostics$upperSpecVarRatio < 8)) 
                    ){
                thinFac <- thinFac/ctrl$thinIncreaseFactor
            }                
            if( thinFac != 1){
                newThin <- as.integer(round(getThin(object) * thinFac))
                if( thinFac > 1){
                    if( newThin > 1.5*ctrl$maxThin) stopDemacConvergenceProblems("too much autocorrelation. Would need to increase thinning interval largely above ",ctrl$maxThin," to ",newThin )
                    newThin <- as.integer(min(newThin, ctrl$maxThin))
                    if( isThinLogs ){
                        thinnedLogs <- squeeze(getSampleLogs(object),fraction=1/thinFac  )  
                        sampleLogs(object) <- thinnedLogs
                        #getNSamplePopulations(thinnedLogs)
                        #plot(asMcmc(thinnedLogs),smooth=FALSE)
                    }
                } 
                thin(object) <- newThin
                message("sampleSann: adjusting thinning interval to ",object@thin)
            }
            object
        })



if(!exists(".isConverged")) setGeneric(".isConverged", function(object,...) standardGeneric(".isConverged"))
#' @export
setMethod(".isConverged", signature="BatchSampler", function(object
                ### report wheter last batch fulfilled convergence criteria
                ,diagnostics=getConvergenceDiagnostics(object)
                ,ctrl=getConvergenceControls(object)
        ){
            isConverged <- if( !object@isBurnin ){
                        # finish early if all chains in Gelman
                        (max(diagnostics$gelmanDiagChains) < ctrl$maxGelmanDiag) &&
                                (diagnostics$gelmanDiagPops < ctrl$maxGelmanDiag) &&
                                !diagnostics$isLogDenDrift 
                    } else { # burnin
                        # finish early only if temperature does not decrease
                        # relaxed convergence criterion for chains
                        nChain <- getNChainInPopulation(getSampleLogs(object))
                        (max(diagnostics$gelmanDiagChains) < sqrt(nChain)*ctrl$maxGelmanDiag) &&
                        (diagnostics$gelmanDiagPops < 1+2*(ctrl$maxGelmanDiag-1)) &&
                        !diagnostics$isLogDenDrift
                    }
            isConverged
        })

if(!exists(".showFinishEarlyMessage")) setGeneric(".showFinishEarlyMessage", function(object,...) standardGeneric(".showFinishEarlyMessage"))
#' @export
setMethod(".showFinishEarlyMessage", signature="BatchSampler", function(object
                ### display a message that sampleBatches finished early
        ){
            diagnostics=getConvergenceDiagnostics(object)
            finishEarlyMsg <- paste("sampleBatches: ",object@cntConsequtiveMinorDecreaseBatches," consequtive batches wich sufficient convergence"
                    ," and no drift in logDensity. Finishing early.",sep="") 
            message( finishEarlyMsg)
            message(paste("gelmanDiagPops=",signif(diagnostics$gelmanDiagPops,2)
                            ," max(gelmanDiagChains)=",signif(max(diagnostics$gelmanDiagChains),2)
                            ," upperSpecVarRatio=",signif(diagnostics$upperSpecVarRatio,2)
                            #," T0-1=",signif(diagnostics$T0-1,2)*100,"%"
                            #," logDen=",paste(signif(lDenLastPart,3),collapse=",")
                            #," T=",paste(signif(TCurr,2),collapse=",")
                            , sep=""))
            
        })

#' @export
computeConvergenceDiagnostics <- function(
        ### compute several diagnostics to check convergence of the chains and populations
        object          ##<< BatchSampler object
        , logsTail      ##<< sample logs to check, usually the tail of the current sample
        , nEffectiveParameter = sqrt(getNParameterWithProposal(logsTail)) ##<< number of effective parameters. Due to correlations it is lower than the number of parameters
){
    nPop <- getNPopulation(logsTail)
    # combining chains, need to slice in order to compare the drift
    logsTailStacked <- stackChainsInPopulation(logsTail, mergeMethod="slice")
    logDenBlocksPops <- lapply( 1:nPop, function(iPop){
                logsPop <- getSampleLogOfPopulation(logsTailStacked, iPop)
                logDenBlocks <- t(adrop( computeBlockLogDensities( logsPop, logsTailStacked )[,,1L ,drop=FALSE],3L))
            })
    if( getNSamplePopulations(logsTail)[1L] < 5L ){
        specVarRatioPops <- rep(NA_real_, nPop)
        gelmanDiagPops <- NA_real_
        gelmanDiagChains <- rep(NA_real_, nPop)
        isLogDenDrift <- NA
        logDenBlocksMean <- sapply(logDenBlocksPops, colMeans)
        logDenBlocksStartEnd <- rbind( logDenBlocksMean, logDenBlocksMean)
    } else {
        specVarRatioPops <- sapply( 1:nPop, function(iPop){
                    parmsPop <- getParametersForPopulation(logsTail, iPop )
                    specVarRatioChains <- apply( parmsPop, 3L, function(parmsChain){
                                #matplot(t(parmsChain), type="l")
                                spec <- spectrum0.ar(t(parmsChain))$spec
                                varSample <- apply(parmsChain, 1, var)
                                spec/varSample
                            })
                    apply(specVarRatioChains, 1, max)    # maximum over chains
                })
        #
        mcl <- asMcmc(logsTailStacked)
        # use try because # cholesky decomposition in gelman.diag may throw errors    
        gelmanDiagPops <- {tmp<- try(gelman.diag(mcl, autoburnin=FALSE), silent=TRUE); if( inherits(tmp,"try-error")) 999 else if(length(tmp$mpsrf)) tmp$mpsrf else tmp$psrf[1] }
        # 
        # convergence in each population
        gelmanDiagChains <- sapply( 1:nPop, function(iPop){
                    logsTailPop <- subsetPopulations(logsTail,iPop)
                    mcl <- asMcmc(logsTailPop)
                    tmp <- try(gelman.diag(mcl, autoburnin=FALSE), silent=TRUE)
                    if( inherits(tmp,"try-error")) 999 else if(length(tmp$mpsrf)) tmp$mpsrf else tmp$psrf[1]
                })
        #
        # drift in resLogDen
        #logDenBlocks <- logDenBlocksPops[[2]]
        maxDrift <- qchisq(object@convergenceControls$qLogDiff, nEffectiveParameter)/2   ##<< difference in LogDensity, below which no drift is signalled
        resLogDenDriftPops <- lapply(logDenBlocksPops, function(logDenBlocks){
                    #matplot(t(logDenBlocks))
                    detectLogDenDrift(logDenBlocks, maxDrift=maxDrift)  
                })
        logDenBlocksStartEnd <- sapply(resLogDenDriftPops, function(resLogDenDrift){
                    attr(resLogDenDrift,"resWTest")[1:2,1L]              
                })
        #resLogDenDrift <-resLogDenDriftPops[[1]] 
        isLogDenDrift <- any(unlist(resLogDenDriftPops))
    }
    #
    diagnostics <- list(
            logDenBlocksPops=logDenBlocksPops      ##<< list of temperated logDensity Components
            ,logDenBlocksStartEnd=logDenBlocksStartEnd  ##<< numeric array (2xnPop) of mean logDenBlocks of first and fouth quartile
            ,gelmanDiagChains=gelmanDiagChains        ##<< numeric vector (nStream) within population gelman diagnostics
            ,gelmanDiagPops=gelmanDiagPops               ##<< gelman diag calculated between populations
            ,specVarRatioPops=specVarRatioPops         ##<< numeric vector (nStream) ratio of spectral density to sample Variance (max over all parameters). This is a measure that increases with autocorrelation.
            ,upperSpecVarRatio=quantile(specVarRatioPops, probs=object@convergenceControls$qSpecVarRatio, na.rm=TRUE)
            #,logDenLastPart=logDenLastPart         ##<< sum of (temperated) block logDensities
            ,isLogDenDrift=isLogDenDrift            ##<< TRUE if logDensity is still relatively and absolutely decreasing between first and fourth quarter of the chain
    )
    diagnostics
}

detectLogDenDrift <- function(
        ### check whether first quartile all the logDensities is significantly smaller than last quartile 
        logDenBlocks	##<< numeric array (nStep x nBlock): logDensity (highest are best)
        , alpha=0.05	##<< the significance level for a difference
        , maxDrift=3*0.5	##<< difference in LogDensity, below which no drift is signalled, should increase with number of parameters
){
    ##details<<
    ## Because of large sample sizes, very small differences may be significantly different.
    ## Use argument \code{maxDrift} to specify below which difference a significant difference is not regarded as drift.
    nr <- nrow(logDenBlocks)
    if( nr < 4L ){
        ds1 <- logDenBlocks[1L,   ,drop=FALSE]        # first quater of the dataset
        ds4 <- logDenBlocks[nr,  ,drop=FALSE]     # last quarter of the dataset
    } else {
        nr4 <- nr %/% 4
        ds1 <- logDenBlocks[1:nr4,   ,drop=FALSE]        # first quater of the dataset
        ds4 <- logDenBlocks[(nr+1-nr4):nr,  ,drop=FALSE]     # last quarter of the dataset
    }
    #iBlock <- 1L
    # matplot(logDenBlocks); matplot(ds1, col="red", add=TRUE); matplot((nr+1-nr4):nr, ds4, col="blue", add=TRUE)
    nBlock <- ncol(logDenBlocks)
    tres <- sapply( 1:nBlock, function(iBlock){
                rs1 <- ds1[,iBlock]
                rs4 <- ds4[,iBlock]
                #resTTest <- try(t.test(rs4,rs1,"greater"), silent=TRUE)
                resWTest <- try(suppressWarnings(wilcox.test(rs4,rs1,"greater", silent=TRUE)))   #likelihood of rs4 is larger than that of rs1
                if( inherits(resWTest,"try-error") ){  # logDen essentially constant
                    #bo <- (diff(resTTest$estimate) >= maxDrift ) && (resTTest$p.value <= alpha)
                    c(mean(rs1),mean(rs4),p=1)
                } else {
                    c( mean(rs1),mean(rs4), p=resWTest$p.value)
                }   
            })
    #tresi <- tres[,1]
    boL <- apply( tres, 2, function(tresi){  
                (diff(tresi[1:2]) >= maxDrift ) && (tresi["p"] <= alpha)
            })
    ret <- any( boL )
    attr( ret, "resWTest") <- tres
    attr( ret, "logDenComp") <- colSums(ds4)       # average logDensity across all blocks of the last quarter
    ##value<< TRUE if any of the logDensities are a significantly greater in the fourth quantile compared to the first quantile of the samples
    ## , attribute \code{resWTest}: numeric maxtrix (logDenStart, logDenEnd, pTTest x nDInfo )
    ## , attribute \code{logDenComp}: numeric vector (nResComp): average logDenT of all components of last quater
    ret
}



#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("BatchSampler")

if(!exists("getThin")) setGeneric("getThin", function(object) standardGeneric("getThin"))
#' @export
setMethod("getThin", signature="BatchSampler", function(object) {object@thin})
if(!exists("thin<-")) setGeneric("thin<-", function(object,value) standardGeneric("thin<-"))
#' @export
setReplaceMethod("thin", signature=c("BatchSampler", "integer"), function(object, value) {object@thin <- value; object})

if(!exists("getConvergenceDiagnostics")) setGeneric("getConvergenceDiagnostics", function(object,...) standardGeneric("getConvergenceDiagnostics"))
#' @export
setMethod("getConvergenceDiagnostics", signature="BatchSampler", function(object,...
        ### Getter method for slot convergenceDiagnostics
        ) {object@convergenceDiagnostics})

if(!exists("getBatchLog")) setGeneric("getBatchLog", function(object,...) standardGeneric("getBatchLog"))
#' @export
setMethod("getBatchLog", signature="BatchSampler", function(object,...
        ### Getter method for slot batchLog
        ) {object@batchLog})

if(!exists("getConvergenceControls")) setGeneric("getConvergenceControls", function(object,...) standardGeneric("getConvergenceControls"))
#' @export
setMethod("getConvergenceControls", signature="BatchSampler", function(object,...
        ### Getter method for slot convergenceControls
        ) {object@convergenceControls})
if(!exists("convergenceControls<-")) setGeneric("convergenceControls<-", function(object,value) standardGeneric("convergenceControls<-"))
#' @export
setReplaceMethod("convergenceControls", signature=c("BatchSampler", "list"), function(object, value
        ### Setter method for slot convergenceControls
        ) {object@convergenceControls <- value; object})


