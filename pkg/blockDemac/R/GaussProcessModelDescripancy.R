# ' @importFrom MCMCpack rinvgamma # ported to this package because of dependency problems (invGamma.R) 
#' @importFrom mvtnorm rmvnorm

#' @export
compileDiscrepancyBlocksForStream <- function(
        ### compile the intermediates and blocks for GP model error of one stream
        x           ##<< locations of observations
        ,obs        ##<< observations
        ,sd2Obs     ##<< variance of observation uncertainty
        , predProcSpec   ##<< named list of one intermediate specification for model prediction
        , predProcSpecEntry=NULL ##<< if predProcSpec returns a list, can provide a name or an index to pick predictions for this data stream 
        , streamSuffix=""	##<< suffix applied to parameters, and specifications to distinguish data streams
        ,priorGammaPsiPars = getGammaParsFromMeanAndVariance(1,0.2) ##<< prior for variance of model discrepancy as a factor of average observation uncertainty
        ,priorInvGammaSd2DiscrPars = getInvGammaParsFromMeanAndMode( 20, 0.05 ) ##<< prior for corellation length
        #,priorInvGammaSd2DiscrPars = c( alpha=0.001, beta=0.001 ) ##<< prior for signal variance
        ,minSd2DiscrFraction = 0 #1/20   ##<< minimum of uncertainty of GP training data as a fraction of model error variance, deprecated: was used to prevent inflating signal variance by  
        ,xPred=numeric(0)  ##<< locations at which to sample model discrepancies, mapped to paramter vector md_1 to md_<nPred>
        ,isSigmaBlockReturned=TRUE  ##<< set to FALSE to omit block for updating sigma_streamSuffix
        #,isUsingExpectedDelta=TRUE  ##<< if TRUE, then deltaO is set to expected value instead of drawn from distribution for performance
        ,isUsingDataBasedSigma=FALSE    ##<< set to TRUE to use an independent Likelihood estimate of signal variance in the calculation of discrepancies, instead of the current signal variance. This reduces performance. 
        ,weightPSd2Discr=length(x)  ##<< the weight for signal variance in the likelihood of correlation length psi
        ,sd2ObsFactor=1             ##<< if sampled discrepancy goes to zero, use this factor to decrease the observation uncertainty during estimation of model discrepancy
        ,nSupportingSet=1L          ##<< integer scalar: number of sets of supporting conditions to examine
        ,isMarkUpdatedPsiBlock=TRUE ##<< if TRUE, block is marked as updated when not accepted, in order to trigger recomputation of supporting locations that involve some randomness
){
    if( !is.list(predProcSpec) || (length(names(predProcSpec)) != 1L) || !is(predProcSpec[[1]],"IntermediateSpecification")) stopDemac(
                "predProcSpec must be a named list with entry of class IntermediateSpecification.")
    if( length(x) != length(obs) ) stopDemac("locations x and observations obs must be of the same length.")
    GPLimits <- .computeGPLimits(x)
    # names of used intermediates and parameters including suffix
    predProcName <- names(predProcSpec)[1]
    iOName <- paste("iO",streamSuffix,sep="_")
    deltaAName <- paste("deltaA",streamSuffix,sep="_")
    logPsiVar <- paste("logPsi",streamSuffix,sep="_")
    logSd2DiscrVar <- if( isSigmaBlockReturned ){
        paste("logSd2Discr",streamSuffix,sep="_")
    } else {
        "logSd2Discr"
    }
    ksiVars <- if(length(xPred)) paste("ksi",seq_along(xPred),streamSuffix,sep="_") else character(0)
    #deltaPredVars <- paste("md",seq_along(xPred),streamSuffix,sep="_")
    meanSd2ObsClosure <- mean(sd2Obs)   # average observation variance
    # training locations based on psi
    xPredVal <- force(xPred)    # force evaluation to assign to closure
    xVal <- force(x)
    iOSpec <- intermediateSpec(
            function(theta, intermediates ){
                psi=exp(theta[logPsiVar])
                if( !is.finite(psi) && is.finite(theta[logPsiVar]) ) stopDemacInvalidChainState("iOSpec: given logPsi yield infinite psi.")
                .computeGPTrainingLocations(psi=psi, x=xVal, minXSpacing=GPLimits$minXSpacing, maxXSpacing=GPLimits$maxXSpacing, xPred=xPredVal 
                        ,nSet = nSupportingSet 
                )
            }
            ,argsUpdateFunction=list()
            ,parameters = c(logPsiVar)
            ,intermediates = character(0) # the default
    )
    deltaASpec <- intermediateSpec(
            function(theta, intermediates, xVal, obsVal, sd2ObsVal, minSd2DiscrFractionVal, priorInvGammaParsVal ){
                predProc <- intermediates[[predProcName]]
                if( length(attr(predProc,"problems")) ) return( structure(list(NA),problems=attr(predProc,"problems")))
                if( length(predProcSpecEntry) ){
                    predProc <- predProc[[predProcSpecEntry]]  
                } 
                if( !length(predProc) ) return(structure( list(NA), problems=paste("deltaASpec: got empty prediction vector")))
                .computeModelDiscrepancies(
                        psi=exp(theta[logPsiVar])
                        , sd2DiscrFac=exp(theta[logSd2DiscrVar])
                        , iORes=intermediates[[iOName]]
                        , predProc=predProc
                        , x=xVal, obs=obsVal, sd2Obs=sd2ObsVal
                        , meanSd2Obs=meanSd2ObsClosure
                        , xPred=xPredVal
                        , minSd2DiscrFraction = minSd2DiscrFractionVal
                        , priorInvGammaPars = priorInvGammaParsVal
                        #,isUsingExpectedDelta = isUsingExpectedDelta
                        , sd2ObsFactor=sd2ObsFactor
                        #, ... 
                )
            }
            ,argsUpdateFunction=list(xVal=x, obsVal=obs, sd2ObsVal=sd2Obs, minSd2DiscrFractionVal = minSd2DiscrFraction
                    , priorInvGammaParsVal=priorInvGammaSd2DiscrPars
            )
            ,parameters = c(logPsiVar, logSd2DiscrVar)
            #,parameters = c(logPsiVar)
            ,intermediates = c( predProcName, iOName)
    )
    intermediateSpecs = c( predProcSpec, list( iOSpec, deltaASpec))
    names(intermediateSpecs) <- c(predProcName, iOName, deltaAName )
    #
    # Metropolis block for logPsi
    force(priorGammaPsiPars) # will be needed in closure
    namesResLogDenPsi <- paste(c("mD","obs","priorPsi"),streamSuffix,sep="_")
    nX <- length(x) 
    # also must depend on a parameter that is updated each cycle, so that 
    # logDen of current state is recalculated based on changed delta
    logDenLogPsiSpec <- blockSpec(logPsiVar, c(logPsiVar),    
            new("MetropolisBlockUpdater",
                    fLogDensity=function(theta, intermediates, logDensityComponents, obsVal, sd2ObsVal){
                        predProc <- intermediates[[predProcName]]
                        if( length(predProcSpecEntry) ) predProc <- predProc[[predProcSpecEntry]]
                        res <- .computeLogDenLogPsi( logPsi=theta[logPsiVar]
                                #, sd2Discr = exp(theta[logSd2DiscrVar])*meanSd2ObsClosure
                                , sd2Discr = intermediates[[deltaAName]]$sd2Discr
                                , predProcX=predProc[1:nX]
                                , deltaAIntermediate=intermediates[[ deltaAName]]
                                , obs=obsVal, sd2Obs=sd2ObsVal
                                , priorGammaPsiPars=priorGammaPsiPars
                                ,weightPSd2Discr=weightPSd2Discr
                        )
                        names(res) <- namesResLogDenPsi # to distinguish from other streams
                        res
                    }
                    ,argsFLogDensity=list(obsVal=obs, sd2ObsVal=sd2Obs)   #priorPsiPars by closure
                    ,logDensityComponentNames = namesResLogDenPsi
                    ,isMarkUpdatedWhenNotAccepted = isMarkUpdatedPsiBlock        # to trigger recomputation of supporting locations
            )
            ,intermediatesUsed=c(predProcName, deltaAName)
    )
    # Gibbs InverseGamma block for sd2Discr
    logSd2DiscrSpec <- blockSpec(logSd2DiscrVar,logSd2DiscrVar,   
            new("FunctionBasedBlockUpdater"
                    , fUpdateBlock=function(...){.updateLogSigmaDiscrByGibbsSamplingGamma(...)}  # allows tracing called function   
                    , argsFUpdateBlock=list( priorInvGammaPars=priorInvGammaSd2DiscrPars, intermediateName=deltaAName
                    , sd2Obs=structure(meanSd2ObsClosure,names=deltaAName)
                    ,isUsingDataBasedSigma=isUsingDataBasedSigma     
                )  
            #, argsFUpdateBlock=list( priorInvGammaPars=priorInvGammaSd2DiscrPars, intermediateName=rep(interMediateNames[3],2), sdObs=rep(structure(meanSd2ObsClosure,names=interMediateNames[3]),2))  
            )        
            ,intermediatesUsed=deltaAName
    )
    # Sampling model discrepancies at other locations
    sampleModelPlusDiscrepancies <- function( 
            ### update Variance by sampling from a scaled inverse Chi-square distribution with using prior information on sigma 
            theta			##<< numeric vector (nParm): current state of parameters used by this block
            ,iParms=seq_along(theta)    ##<< integer vector (nParmToUpdate): index of parameters to update within theta
            ,upperParBounds ##<< named numeric vector, upper limits of the parameters to update  
            ,lowerParBounds ##<< named numeric vector, lower limits of the parameters to update  
            ,intermediates=list() ##<< intermediate result for current state xC, see end of vignette on using two Densities
    # by closure,IGPars0            ##<< vector with entries "alpha" and "beta"
    ){
        #print("recover at sampleModelDiscrepancies"); recover()
        deltaA <- intermediates[[ deltaAName ]]  # actual sampling of deltaP and ksiP was done together with deltaO in intermediate
        ##value<< list with components
        list(	##describe<<
                isUpdated=TRUE		##<< boolean, if changed block parameters, always true with Gibbs sampling
                , xC=deltaA$ksiP 	##<< numeric vector: components of position in parameter space that are being updated
        )	##end<<
    }
    ksiSpec <- blockSpec(ksiVars,c(ksiVars),   
            new("FunctionBasedBlockUpdater"
                    , fUpdateBlock=sampleModelPlusDiscrepancies
            #, argsFUpdateBlock=list( IGPars0=priorInvGammaSd2DiscrPars)  
            )        
            ,intermediatesUsed=deltaAName
    )
    #
    blockSpecs <- list(logDenLogPsi=logDenLogPsiSpec, logSd2Discr=logSd2DiscrSpec, ksi=ksiSpec)
    if( length(xPred)==0 ) blockSpecs$ksi <- NULL
    if( !isTRUE(isSigmaBlockReturned) ) blockSpecs$logSd2Discr <- NULL
    names(blockSpecs) <- paste(names(blockSpecs),streamSuffix,sep="_")
    #
    fConstrainingSubSpace <- .getFConstrainingSubSpace(GPLimits$logPsiMin,GPLimits$logPsiMax, logPsiVar )
    ##value<< list with entries
    list(
            intermediateSpecs=intermediateSpecs ##<< list with intermediateSpecifications predProc_suffix, iO_suffix, and deltaA_suffix  
            ,blockSpecs=blockSpecs              ##<< update block specifications logDenLogPsi_suffix, and logSd2Discr_suffix 
            ,GPLimits=GPLimits                  ##<< results of \code{\link{computeGPLimits}}
            ,fConstrainingSubSpace=fConstrainingSubSpace
            ,namesProcessSamples=ksiVars
            ,meanSd2Obs=structure( meanSd2ObsClosure,names=deltaAName) ##<< named numeric scalar: average observation uncertainty, with name of deltaA intermediate for this stream 
    )
}
attr(compileDiscrepancyBlocksForStream,"ex") <- function(){
    if(FALSE){
        require(twMisc) # dumpOnError
        
        #---- defining the model and observations
        set.seed(0815)      # for reproducible results
        simpleModel <- function(theta,xval){ theta[1]*xval }
        sdObsTrue <- 0.01
        thetaTrue <- c(theta=0.65) # the true parameter vector
        nObs <- 11
        #nObs <- 61
        x <- seq(0,4,length.out=nObs+1)[-1]
        obsTrue <- simpleModel(thetaTrue,x)/(1+x/20)
        obs <- obsTrue + rnorm(length(x), sd=sdObsTrue)
        sd2Obs <- sdObsTrue^2
        thetaHat <- structure(coef(lm(obs ~x-1)), names="theta")
        theta0 = c(theta=as.numeric(thetaHat), logPsi=log(0.6), logSd2Discr=log(0.5))
        covarTheta0 = diag(c(theta0*0.2)^2, nrow=length(theta0))    # spread for initial population
        # ---- define the prediction intermediate 
        xPred <- c(1.5,3.8,6)
        nPred <- length(xPred)
        predProcSpec <- intermediateSpec(
                dumpOnError(function(theta, intermediates, xVal ){
                            list(pred = simpleModel(theta["theta"],xVal))
                        })
                ,argsUpdateFunction=list(xVal=c(x,xPred))
                ,parameters = c("theta")
                ,intermediates = character(0)  
        )
        streamSuffix <- "obs1"
        theta0N <- c(theta0, structure(rep(0, length(xPred)),names=paste("ksi",seq_along(xPred),sep="_")))
        names(theta0N)[2:(3+nPred)] <- paste(names(theta0N)[2:(3+nPred)],streamSuffix,sep="_")
        covarTheta0 = diag(c(theta0N[1:3]*0.2,rep(0.2,nPred))^2, nrow=length(theta0N))    # spread for initial population
        #
        res <- compileDiscrepancyBlocksForStream(x, obs, sd2Obs
                        , predProcSpec = list(predSpeed=predProcSpec), predProcSpecEntry = 1L
                        , streamSuffix, xPred=xPred)
        #
        predProcRes0 <- getUpdateFunction(res$intermediateSpecs$predSpeed)(theta0, list(), xVal=c(x,xPred)) 
        iORes0 <- getUpdateFunction(res$intermediateSpecs$iO_obs1)(theta0N, list()) 
        deltaARes0 <- getUpdateFunction(res$intermediateSpecs$deltaA_obs1)(theta0N, list(iO_obs1=iORes0, predSpeed=predProcRes0)
                , xVal=x, obsVal=obs, sd2ObsVal=sd2Obs, minSd2Discr=1/20, priorInvGammaParsVal=getInvGammaParsFromMeanAndMode( 20, 0.05 ))
        #predicting model discrepancy at other x
        #KKp <- outer(xPred, deltaARes0$xO, .corG, psi=exp(theta0N["logPsi_obs1"]))
        #as.vector(KKp %*% deltaARes0$KyInvY)
        intermediates0 <- list(predSpeed=predProcRes0, iO_obs1=iORes0, deltaA_obs1=deltaARes0)
        testLogDenLogPsi <- getFLogDensity(getBlockUpdater(res$blockSpec$logDenLogPsi_obs1))(
                theta0N
                , intermediates0, numeric(1)
                , obs=obs, sd2Obs=sd2Obs
        #, priorPsiPars=priorPsiPars 
        )
        argsFUpdate <- getArgsFUpdateBlock(getBlockUpdater(res$blockSpec$logSd2Discr_obs1))
        testUpdateSigma <- do.call(getFUpdateBlock(getBlockUpdater(res$blockSpec$logSd2Discr_obs1)), c(list(  
                                theta0N, iParms=match("logSd2Discr_obs1",names(theta0N))
                                , upperParBounds=c(), lowerParBounds=c()
                                , intermediates=intermediates0
                        ),argsFUpdate) ) 
        testSampleKsi <- getFUpdateBlock(getBlockUpdater(res$blockSpec$ksi_obs1))(  
                theta0N, iParms=grep("^ksi_.*_obs1$",names(theta0N))
                , upperParBounds=c(), lowerParBounds=c()
                , intermediates=intermediates0
        )
        #
        # Metropolis block sampling theta, based on intermediate predProcSpec
        logDenTheta <- function( theta, intermediates, logDensityComponents, obs, sd2Obs ){
            predProc <- intermediates$predSpeed[[1]] 
            deltaA <- intermediates$deltaA_obs1$deltaA # heed the streamSuffix
            logPriorInvGammaSd2Discr <- intermediates$deltaA_obs1$logPriorInvGammaSd2Discr # heed the streamSuffix
            residuals <- obs - (predProc[seq_along(x)] + deltaA) 
            #sd2Obs <- exp(theta["logSd2Obs"])
            return( c(obs1= -1/2 * sum( ((residuals)^2)/sd2Obs ), md_obs1=logPriorInvGammaSd2Discr) )
        }
        testLogDenTheta <- logDenTheta(theta0, intermediates0, numeric(1), obs=obs, sd2Obs=sd2Obs )
        logDenThetaSpec <- blockSpec("theta",c("theta","logPsi_obs1","logSd2Discr_obs1"),    
                new("MetropolisBlockUpdater",
                        fLogDensity=function(...){logDenTheta(...)},
                        argsFLogDensity=list( obs=obs, sd2Obs=sd2Obs),
                        logDensityComponentNames = names(testLogDenTheta)
                )
                ,intermediatesUsed=c("predSpeed","deltaA_obs1")
        )
        #
        res <- compileDiscrepancyBlocksForStream(x, obs, sd2Obs
                , predProcSpec = list(predSpeed=predProcSpec), predProcSpecEntry = 1L
                , streamSuffix=streamSuffix, xPred=xPred
                ,minSd2DiscrFraction = 0 #1/20
        )
        intermediateSpecs <- c(res$intermediateSpecs)
        blockSpecs <- c(list(bTheta=logDenThetaSpec), res$blockSpecs)
        #mtrace(logDenTheta, recover)   #untrace(logDenTheta)
        sampler <- newPopulationSampler( blockSpecs, theta0N, covarTheta0
                ,intermediateSpecifications=intermediateSpecs
                , subSpace=res$fConstrainingSubSpace(new("SubSpace"))
        )
        
        sampler <- setupAndSample(sampler, nSample=60L)
        sampler <- setupAndSample(sampler, nSample=60L)
        #sampler <- setupAndSample(sampler, nSample=8L)
        #sampler <- setupAndSample(sampler, nSample=200L)
        sampleLogs <- getSampleLogs(sampler)
        plot( asMcmc(sampleLogs), smooth=FALSE )
        stacked <- adrop(getParametersForPopulation(stackChainsInPopulation(subsetTail(sampleLogs, 0.6)),1L),3)
        summary(t(stacked))
        plot(density(stacked["theta",]))
        
        thetaMean <- rowMeans(stacked)
        thetaCf <- apply(stacked, 1, quantile, probs=c(0.05,0.95))
        yPredTrue <- simpleModel(thetaTrue,xPred)/(1+xPred/20)
        yPredSampled <- thetaMean[3+(1:length(xPred))]
        cfYPredSampled <- thetaCf[,3+(1:length(xPred)), drop=FALSE]
        plot(obsTrue ~ x, col="blue", xlim=range(c(x,xPred)), ylim=range(c(obs,yPredTrue,yPredSampled)))
        points(obs ~ x, col="green")
        #lines(simpleModel(thetaHat,x) ~ x)
        lines(simpleModel(thetaMean["theta"],x) ~ x)
        points(yPredTrue  ~ xPred, col="blue")
        points(yPredSampled ~ xPred, col="orange")    
        arrows(xPred,cfYPredSampled[1,], xPred, cfYPredSampled[2,], angle=90, code=3, length=0.05, col="orange")
        legend("topleft",inset=c(0.01,0.01)
                , legend=c("observations","true values","model prediction","process prediction")
                ,col=c("green","blue","black","orange"), pch=c("o","o","","o")
                ,lty=c(0,0,1,1)
        )
    }
}


.computeGPLimits <- function(
        ### computing limits on estimates of correlation length psi and spacing of trainign points for Gaussian Process
        x                   ##<< locations of observations
        , maxKernelDim=60   ##<< maximum of dimensionality of the Kernal matrix
){
    ##details<< 
    ## For large psi, some matrices become indeterminate, hence bound by range of x-values.
    ## For very small psi, model error becomes zero between supporting points. However, 
    ## we want to have a smooth model error between supporting points, else we model the noise.
    ## Hence, bound psi by only a bit lower than minimal allowed average distance between supporting points.
    ## If every observation becomes a supporting point, then we model the noise.
    ## Hence, we ensure not allow smaller distance between supporting points than 4 times the average distance
    ## covered by the observations.
    xRange <- diff(range(x))
    maxXSpacing <- xRange / (3+2) 
    minXSpacing = max( xRange/maxKernelDim, min( maxXSpacing, xRange/length(x)*4))
    #logPsiMin = log(xRange/(length(x)/2)) 
    logPsiMin = log(minXSpacing / 1.5)   # not much smaller than minimal xSpacing
    logPsiMax = log(xRange)   # for larger psi, inversion of matrix becomes singlular
    ##value<< list with entries
    list(
        logPsiMin = logPsiMin            ##<< minimum correlation length about halv the distance between supporting points
        ,logPsiMax = logPsiMax           ##<< maximum correlation length of about the range of x
        ,minXSpacing = minXSpacing       ##<< minimum distrance of supporting points so that there are about 2 to 3 points in between
        ,maxXSpacing = maxXSpacing              ##<< maximum distance of training points so that there are at least 3 points in addition to start and end 
        ,psi0 = min(exp(logPsiMax),5*exp(logPsiMin))    ##<< initial estimate of correlation length, of about 10 oints
    )
}

# need to export to work on cluster
#' @export
.corG <- function(x1,x2,psi){
    exp(-((x1-x2)/psi)^2)
}

# need to export to work on cluster
#' @export
.computeGPTrainingLocations <- function(
        ### compute subsets of locations, xO, where observation error will be sampled from Gaussian Process
        psi                 ##<< current estimate of correlation length of Gaussian Process
        ,x                  ##<< locations, e.g. time, at witch observations are available
        , minXSpacing       ##<< minium distance between neighboring xO
        , maxXSpacing       ##<< maximum distance between neighboring xO
        , xPred=numeric(0)  ##<< numeric vector: locations where predictions of GP are recorded
        , psiSpacingFac=1.5 ##<< multiple of psi, of xSpacing
        , avgFrac=1/10      ##<< length of the range over which to average observations, as fraction of distance between supporting locations
        , nSet=1L           ##<< number of alternative sets of supporting locations returned
){
    # determine subset of xO by taking points all every ~psi locations
    if( !is.finite(psi) ) stopDemac("provided a psi that is not finite. Check variable nameing in theta.")
    xOSpacing <- xOSpacing0 <- min(max((psiSpacingFac*psi), minXSpacing), maxXSpacing)   # to prevent too large matrices do not go beyond a minimum of spacing between x
    xd0 <- xd <- c(0,cumsum(diff(x)))   # distances to x[1]
    # choose the second point o2 among those locations that are less than 1/4 away from x0Spacing
    # keep the first points, else GP goes to zero at the edge
    o2Min <- which.min(abs(xd0 - xOSpacing*3/4))
    o2Max <- which.min(abs(xd0 - xOSpacing*5/4))
    nChooseFrom <- (o2Max-o2Min)+1
    o2Offsets <- sample.int(nChooseFrom, size=nSet, replace=TRUE)-1
    # introduce some jitter to select different locations by varying the distance between supporting locations
    xOSpacings <- matrix( runif( nSet* (2L+ceiling(diff(range(x))/xOSpacing)), min=xOSpacing*0.9, max=xOSpacing*1.1 ), ncol=nSet )
    xOTargetCums <- apply(xOSpacings, 2, cumsum)
    #iSet <- 1L
    ans <- lapply( 1:nSet, function(iSet){
        o2 <- o2Min + o2Offsets[iSet]
        xd <- xd0[-(1:o2)]-xd0[o2]  # distances of other points to o2
        xOTargetCum <- xOTargetCums[,iSet]
        i <- 0L; j <- 1L
        iO2 <- integer(length(xd)) 
        while( j < length(xd)){
            i <- i+1L
            jn <- j+which.min( (xd[-(1:j)]-xOTargetCum[i])^2 )
            if( !length(jn) ) stop(".computeGPTrainingLocations: encontered jn of zero length.") 
            j <- jn
            iO2[i] <- j
        }
        iO <- c(1L, o2, o2+iO2[1:i])
        # make sure last point is in set
        #nX <- length(x)
        #mIO <- max(iO)
        #if( x[nX] > (x[mIO]+xOSpacing/2) ) iO <- c(iO,nX) else iO[length(iO)] <- nX
        iAvg <- rep(NA_real_, length(x))
        for( i in seq_along(iO)){
            ii <-  ( if(i==1) TRUE else (x > x[iO[i-1L]]) ) &    # only consider x above previous supporting location
                   ( if(i==length(iO)) TRUE else (x < x[iO[i+1L]]) ) &  # below next supporting location
                   abs(x - x[iO[i]] ) < xOSpacing*avgFrac/2 
            iAvg[ii] <- i    
        }
        if( length(na.omit(unique(iAvg))) != length(iO) ){ print("encountered problems averaging locations for supporing points."); recover() }
        xO <- x[iO]
        #plot(1*x ~ x, pch="+", col="grey"); points(1*xO~xO, col="red")
        xOP <- c(xO,xPred)
        LambdaOO <- outer(xO, xO, .corG, psi=psi) 
        LambdaRO <- outer(x[-iO], xO, .corG, psi=psi)
        #LambdaRR <- outer(x[-iO], x[-iO], .corG, psi=psi)
        #LambdaPP <- outer(xOP, xOP, .corG, psi=psi)
        #LambdaPO <- outer(xOP, xO, .corG, psi=psi)
        LambdaPP <- outer(xPred, xPred, .corG, psi=psi)
        LambdaPO <- outer(xPred, xO, .corG, psi=psi)
        ##value<< list of nSet sets of supporting locaitons. Each is a lists with entries     
        list(
                iO = iO                 ##<< indices within x, for which model discrepancies are sampled
                ,iAvg = iAvg            ##<< index vector (length(x)): giving indices of locations to average observations for training the GP
                ,LambdaOO = LambdaOO
                ,LambdaRO = LambdaRO
                #,LambdaRR = LambdaRR
                ,LambdaPP = LambdaPP    ##<< correlations among c(x,xPred)
                ,LambdaPO = LambdaPO
        )
    })
}

.tmp.f <- function(){
    plot( seq_along(x) ~ x )
    abline(v=x[c(1,o2)])
    abline(v=x[iO], col="green")
    
    xR <- x[-(1:o2)]-x[o2]
    plot(xd ~ xR )
    abline(v=xR[iOR])
    abline(v=xR[iOR], col="red")
    

}

# need to export to work on cluster
#' @export
.getFConstrainingSubSpace <- function(
        ### provide a function that constrains a subspace for given logPsiMin and logPsiMax
        logPsiMin           ##<< numeric scalar: lower bound 
        ,logPsiMax          ##<< numeric scalar: upper bound
        ,varName="logPsi"	##<< name of the parameter to constrain
){
    ##value<< a function that constrains a subspace to the given bounds
    constrainSubSpace <- function(
            ### constrain subspace low lower and upper parameter bound
            subSpace=new("SubSpace")
    ){
        getSplitSubSpaces(
            getSplitSubSpaces(subSpace, structure(logPsiMax,names=varName) )$lower
            ,structure(logPsiMin,names=varName)
        )$upper
    }
}
attr(.getFConstrainingSubSpace,"ex") <- function(){
    fConstr <- .getFConstrainingSubSpace(0.5,1.5,"logPsi")
    cSubSpace <- fConstr(new("SubSpace"))
    # add this subspace to the sampler: by providing it to newSampler 
}


# need to export to work on cluster
#' @export
.computeModelDiscrepancies <- function(
        ### compute model discrepancies MD given locations, sigma, psi and model discrepancies
        #sd2DiscrFac    ##<< estimated variance of factor of MD/sdObs 
        psi        ##<< estimated correlation length
        ,sd2DiscrFac    ##<< numeric scalar: current signal variance (normalized by meanSd2Obs) 
        ,iORes          ##<< indices within x to sample from gaussion process GP
        ,predProc   ##<< intermediate giving model predictions at c(x, xPred)
        ,x          ##<< locations, e.g. time, of observations 
        ,obs        ##<< observations
        ,sd2Obs     ##<< numeric scalar: variance of observation uncertainty
        ,meanSd2Obs      ##<< numeric scalar: mean of observation uncertainty variance
        ,minSd2DiscrFraction = 1/20   ##<< minimum of uncertainty of GP training data as a fraction of model error variance
        ,xPred=numeric(0)  ##<< locations at which to sample model discrepancies, mapped to paramter vector md_1 to md_<nPred>
        ,priorInvGammaPars ##<< prior of signal variance s2Discr #= getInvGammaParsFromMeanAndMode( 0.3^2, 0.2^2 )
        ,sd2ObsFactor=1       ##<< if sampled discrepancy goes to zero, decrease the estimate of observation uncertainty by lowering this factor
){
    #isUsingMLESigma=TRUE  ##<< set to TRUE, if estimate sd2Discr by maximising likelihood instead of expecting misifit of n    
    if( length(sd2Obs) == 1L ) sd2Obs <- rep(sd2Obs, length(obs))
    #if( !is.finite(sd2DiscrFac) ) stopDemac("provided a sd2Discr that is not finite. Check variable nameing in theta.")
    nX <- length(x)
    nP <- length(xPred)
    z <- (obs - predProc[1:nX])     # actual differences
    
    #iOResI <- iORes[[1L]]
    resSets <- lapply( iORes, function(iOResI){
        iO <- iOResI$iO
        xO <- x[iO]
        nO <- length(iO)
        iAvg <- iOResI$iAvg
        #
        ##details<< 
        ## In order to diminish dependence of the discrepancy on the 
        ## choice of supporting locations and their random observation error,
        ## Observations and their variance are over neighboring locations. 
        zO <- tapply(z,iAvg, mean)
        sd2ObsO <- as.vector(tapply(sd2Obs, iAvg, function(sd2){sum(sd2)/length(sd2)^2}))
        sd2Discr <- sd2DiscrSampled <- as.vector(sd2DiscrFac * meanSd2Obs)
        #
        KOO <- sd2Discr * iOResI$LambdaOO
        KRO <- sd2Discr * iOResI$LambdaRO
        ##details<<
        ## If the model error is fluctuating, the discrepancies and their variance
        ## may be estimated too low and go to zero.
        ## A workaround is using a smaller observation variance when computing expected discrepancies
        ## (argument sd2ObsFactor)
        Ky <- KOO + diag(sd2ObsO*sd2ObsFactor, nrow=length(xO))
        if( nP != 0 ){
            xOP <- c(xO,xPred)  # prediction at xO and xP
            if( length(xOP) > length(predProc)) stopDemac(
                   "predProc intermediate does not provide enough model predictions (maybe does not including xPred?)")
            Kpp <- sd2Discr * iOResI$LambdaPP
            Kpo <- sd2Discr * iOResI$LambdaPO
            KySolved <- solve(Ky, cbind(zO, t(Kpo) )) # call solve only once
            solvedKyZ <- KySolved[,1]
            solvedKyKop <- KySolved[,1+seq_along(xPred)]
            #
            ##details<<
            ## When predicting model discrepancies at prediction locations,
            ## They are sampled from the GP. Also if the discrepancy at
            ## supporting locations is set to the expected value
            deltaPMu <- as.vector(Kpo %*% solvedKyZ)  # expected value at prediction locations
            deltaPCov <- Kpp - Kpo %*% solvedKyKop
            deltaPCovSymm <- (deltaPCov + t(deltaPCov)) / 2     # else sometimes rmvnorm does complain 
            deltaP <- pred <- 
                    as.vector(rmvnorm(1, deltaPMu, deltaPCovSymm, method="svd")) # realization
            ksiP <- predProc[nX+(1:nP)] + deltaP
            if( any(!is.finite(ksiP))) stopDemac("encountered non-finite process predictions")
        } else {
            solvedKyZ <- solve(Ky, zO)
            deltaP <- ksiP <- numeric(0)
        }
        deltaA <- numeric(length(x))        # observation errors of all x
        deltaO <- deltaA[iO] <- deltaOMu <- as.vector(KOO %*% solvedKyZ)  # maybe smoothed compared to zO
        # in computing deltaR, do not account for observation uncertainty (use KOO and instead of Ky), because using deltaO instead of y 
        # but if deltaO is the expected value, then  solvedKooDeltaO corresponds to solvedKyZ
        solvedKooDeltaO <- solvedKyZ
        deltaR <- deltaA[-iO] <- as.vector(KRO %*% solvedKooDeltaO)
        SSQpred <- sum( (z - deltaA)^2/sd2Obs )
        #------- plotting discrepancies and model
        # plot(deltaA ~ x, type="l", ylim=range(c(z)) );    points( z ~ x, pch="+", col="grey");  points(z[iO] ~ x[iO], col="red", pch="+"); points(zO ~ x[iO], col="red"); abline(h=0, col="grey");    points(deltaOMu ~ xO, col="maroon" );    points(deltaO~xO,pch="o", col="green");  
        # plot(deltaA ~ x, type="l" );   points(deltaOMu ~ xO, col="maroon" );    points(deltaO~xO,pch="o", col="green")   
    #if( nX > 20 ) recover()
        ##details<<
        ## The product delta %*% K^-1 %*% delta is approximated by delta_S %*% KSS^-1 delta_S.
        ## Otherwise the inversion of the K matrix involves an inversion of a submatrix of dimension n_R
        ## Additionally, this becomes numerically singular easily.
        ## With some preliminary testing (both with sparse and rich streams), the approximation was precise 
        ## within a tolerance of 1permill of the original product.
        #solvedKooDeltaO <- solve(KOO, deltaO)
        normDelta <- tmp <- as.vector(deltaO %*% solvedKooDeltaO )      # rely on the approximation by deltaO for all delta
        #invCovDeltaBlocks <- computeInverseOfBlockedMatrix( KOO, t(KRO), KRO, KRR )
        #normDelta <- blockedMatrixNorm( deltaO, deltaR, invCovDeltaBlocks )
        #if( abs(normDelta - tmp) > normDelta/100 )  recover()
        logPDeltaO <- -1/2 * normDelta  
        #
        ##value<< a list with values
        res <- list(
                deltaA=deltaA           ##<< sample of model discrepancies at all observations locations (xO and xR) 
                ,deltaP=deltaP          ##<< sample of model discrepancies at xPred
                ,ksiP=ksiP              ##<< current model predictions + deltaP, i.e. sample of process value
                #,Ky=Ky                  ##<< GP kernal with uncertain observations at all xO (subset of x)
                ,xO=xO                  ##<< the positions of the training data
                ,iO=iO                  ##<< index of subset of training points xO among x
                ,psi=psi                ##<< hyperparameter used in the current sample
                ,sd2Discr=sd2Discr      ##<< hyperparameter used in the current sample
                ,meanSd2Obs=meanSd2Obs  ##<< mean observation variance, used to normalize sd2Discr
                ,zO=zO                  ##<< model-data residual at supporting locations (aggregate of neighboring locations)
                ,sd2ObsO=sd2ObsO        ##<< observation variance of zO (aggregate of neighboring locations)
                #,logPriorInvGammaSd2Discr=logPriorInvGammaSd2Discr  ##<< unnormalized prior probability of this expected signal variance
                ,LambdaOO = iOResI$LambdaOO  ##<< correlation matrix used for deltaO
                ,solvedLambdaOODeltaO = sd2Discr * solvedKooDeltaO ##<< Lambda^-1 deltaO
                #,solvedKyy = KySolved[,1]   ##<< Ky^-1 (obs - pred)
                #,solvedKooDeltaO = solvedKooDeltaO  ##<< KOO^-1 deltaO
                #,KOO = KOO
                ,logPDeltaO = logPDeltaO    ##<< penalty term for model discrepancies in log-Density
                ,SSQpred = SSQpred      ##<< sum( (obs-(pred+delta))^2/sdObs^2 )
        )
    })
    #logPDeltaOs <- sapply( resSets, "[[", "logPDeltaO")
    #quant80 <- sort(logPDeltaOs)[ round(0.8*length(logPDeltaOs)) ]
    #i <- which(logPDeltaOs == quant80)[1]
    ##details<< 
    ## Model discrepancies are computed for several sets of supporting locations.
    ## For each the weighted Sum of Squares of misfits is computed
    ## From all the supporting locations sets, the result is chosen that corresponds to the 80percentile of misfits
    SSQs <- sapply( resSets, "[[", "SSQpred")
    quant80 <- sort(SSQs)[ round(0.8*length(SSQs)) ]
    i <- which(SSQs == quant80)[1]
    ans <- resSets[[i]]
    # the following for testing update of sd2DiscrFac
    #intermediates <- list(deltaA_obs1=res); sd2Obs <- list(deltaA_obs1=meanSd2Obs)
    ans
}

.tmp.f <- function(){
    expectedSd2DiscrFac <- seq(0.01,0.7,by=0.01)
    plot( dinvgamma(expectedSd2DiscrFac, priorInvGammaPars["alpha"], priorInvGammaPars["beta"]) ~ expectedSd2DiscrFac )
    plot( log(dinvgamma(expectedSd2DiscrFac, priorInvGammaPars["alpha"], priorInvGammaPars["beta"])) ~ expectedSd2DiscrFac )
    #
    logPriorInvGammaSd2Discr <- as.vector(log(expectedSd2DiscrFac)*(-priorInvGammaPars["alpha"]-1) - priorInvGammaPars["beta"]/expectedSd2DiscrFac)
    plot(logPriorInvGammaSd2Discr ~ expectedSd2DiscrFac)
    
}

.estimateSignalVariance <- function(
        ### estimate signal variance for given correlations and data
        zO      ##<< numeric vector: model data misfit 
        , LambdaOO=outer(xO, xO, .corG, psi=psi) ##<< correlation matrix
        , sd2ObsO ##<< observation variance (vector of length zO)
        , psi   ##<< alternatively to specifying LambdaOO, correlation length psi can be specified together with xO or spacing ds
        , xO=.computeGPTrainingLocations(ds, psi, x)    ##<< supporting locations of the GP  
        , ds    ##<< distance between supporting locations
        , x     ##<< locations of observations
){
    #solvedLambdaOOZO <- solve(LambdaOO, zO)
    #sd2Discr <- as.vector(zO %*% solvedLambdaOOZO) / length(zO)
    nO <- length(zO)
    invLambdaOO <- solve(LambdaOO)
    sd2Discr1 <- sd2Discr0 <- sd2Discr <- as.vector(zO %*% invLambdaOO %*% zO) / nO
    if( sd2Discr < 0 ){
        stop("encountered negative variance.")
        sd2Discr <- sd2Discr0
    }
    sd2Discr <- if( sd2Discr0 > 6*mean(sd2ObsO)){
        sd2Discr0        
    } else {
        # Correct by part of variance explained by sd2_eps
        sd2DiscrI <- numeric(10)
        sd2DiscrI[1L] <- mean(sd2ObsO)/2
        sd2DiscrI[2L] <- sd2Discr0
        tol=sd2Discr0/20
        sdObsInv <- diag(1/sd2ObsO, nrow=nO)
        i <- 3L
        while( (i <= 10) && (abs(sd2DiscrI[i-1L] - sd2DiscrI[i-2L]) > tol) ){
            sd2Discr <-  (1*sd2DiscrI[i-1L] + 1*sd2DiscrI[i-2L])/2          
            W <- invLambdaOO/sd2Discr + sdObsInv
            invW <- solve(W)
            W2 <- invLambdaOO %*% invW %*% invLambdaOO /sd2Discr
            sd2DiscrI[i] <- (sd2Discr0 -  as.vector(zO %*% W2 %*% zO) / nO)
            i <- i +1L
            
        }
        #plot( sd2DiscrI )
        sd2Discr <-  (1*sd2DiscrI[i-1L] + 1*sd2DiscrI[i-2L])/2
    }
    #Ky <- sd2Discr*LambdaOO + diag(sd2ObsO, nrow=nO); zO %*% solve(Ky,zO)
    sd2Discr
}


.estimateSignalVarianceOpt <- function(
        ### estimate signal variance for given correlations and data
        zO      ##<< numeric vector: model data misfit 
        , LambdaOO=outer(xO, xO, .corG, psi=psi) ##<< correlation matrix
        , sd2ObsO ##<< observation variance (vector of length zO)
        , psi   ##<< alternatively to specifying LambdaOO, correlation length psi can be specified together with xO or spacing ds
        , xO=.computeGPTrainingLocations(ds, psi, x)    ##<< supporting locations of the GP  
        , ds    ##<< distance between supporting locations
        , x     ##<< locations of observations
        , solvedLamdaOOzO = solve( LambdaOO, zO)    ##<< may provide for performance reasons
        #, sd2Discr0 = as.vector(zO %*% solvedLamdaOOzO) / length(zO)    ##<< initial estimate of discrepancy variance
){
    #solvedLambdaOOZO <- solve(LambdaOO, zO)
    #sd2Discr <- as.vector(zO %*% solvedLambdaOOZO) / length(zO)
    nO <- length(zO)
    #x <- seq(2, 2*sd2Discr0, length.out=31)
    #logPZO <- sapply(x, .logPZO, LambdaOO=LambdaOO, sd2ObsO=sd2ObsO, zO=zO  )
    #plot( logPZO ~ x)
    #abline(v=sd2Discr)
    sd2DiscrUpper <- as.vector(zO %*% solvedLamdaOOzO) / length(zO) 
    resOpt <- optimize(.logPZO, interval=c(0, sd2DiscrUpper), maximum = TRUE
            ,tol=sd2DiscrUpper/20
            ,LambdaOO=LambdaOO, sd2ObsO=sd2ObsO, zO=zO
    )
    sd2Discr <- resOpt$maximum
    sd2Discr
}

.logPZO <- function(
        ### Likelihood of zO from Gaussian process 
        sd2Discr
        , LambdaOO  ##<< correlation matrix
        , sd2ObsO   ##<< observation variance (vector of length zO)
        , zO        ##<< observed model-data residuals
){
    nO <- length(zO)
    Ky <- sd2Discr*LambdaOO + diag(sd2ObsO, nrow=nO)
    U <- chol(Ky)
    alpha <- backsolve( t(U), backsolve(U,zO))
    logPZO <- -1/2*zO%*%alpha -sum(log(diag(U))) #-nO/2*log(2*pi) 
}




#' @export
.updateLogSigmaDiscrByGibbsSamplingGamma <- function( 
        ### update Variance by sampling from a scaled inverse Chi-square distribution with using prior information on sigma 
        theta			##<< numeric vector (nParm): current state of parameters used by this block
        ,iParms=seq_along(theta)    ##<< integer vector (nParmToUpdate): index of parameters to update within theta
        ,upperParBounds ##<< named numeric vector, upper limits of the parameters to update  
        ,lowerParBounds ##<< named numeric vector, lower limits of the parameters to update  
        ,intermediates=list() ##<< intermediate result for current state xC, see end of vignette on using two Densities
        ,intermediateNames="deltaA"	##<< name of the intermediate, usually combined with a suffix for given stream
        ,sd2Obs              ##<< named numeric vector: average observation uncertainty. Names must correspond to argument intermediateNames 
        ,priorInvGammaPars #= getInvGammaParsFromMeanAndMode( 0.3^2, 0.2^2 )
        ,isUsingDataBasedSigma=FALSE    ##<< set to TRUE to use an independent Likelihood estimate of signal variance in the calculation of discrepancies, instead of the current signal variance. This reduces performance. 
){
    #intermediateName <- intermediateNames[2]
    #intermediateName <- "deltaA_obs1"
    # deltaO%*%solve(LambdaOO,deltaO) / length(deltaO)      # estimate of signal variance
    # allow using one discrepancy as the mean over estimates of several blocks, denoted by intermediateNames
    sd2DiscrFacStreams <- sapply(intermediateNames, function(intermediateName){
                int <- intermediates[[ intermediateName ]]
                if( length(attr(int,"problems")) ) return(exp(theta))   # do not resample, because have no deltaO
                iO <- int$iO
                solvedLambdaOODeltaO <- int$solvedLambdaOODeltaO
                deltaO <- int$deltaA[iO]
                #LambdaOO <- int$LambdaOO
                #normDeltaO <- tmp <- deltaO%*%solve(LambdaOO,deltaO)
                normDeltaO <- tmp <- deltaO %*% solvedLambdaOODeltaO
                #sd2Discr0 = exp(theta)*sd2Obs[[intermediateName]]   # value used for current delta
                if( isUsingDataBasedSigma ){
                    sd2DiscrData <- .estimateSignalVarianceOpt( int$zO, int$LambdaOO, int$sd2ObsO)
                    KOO <- sd2DiscrData*int$LambdaOO 
                    Ky <- KOO + diag(int$sd2ObsO, nrow=length(iO))
                    solvedKyZ <- solve(Ky, int$zO)
                    deltaOD <- as.vector(KOO %*% solvedKyZ)
                    solvedLambdaDeltaO <- sd2DiscrData * solvedKyZ  # only if deltaOD are expected values
                    normDeltaO <- deltaOD %*% solvedLambdaDeltaO
                }
                #
                IGParsAdd <- c( length(deltaO)/2, 1/2/sd2Obs[[intermediateName]]  * normDeltaO)
                IGPars <- priorInvGammaPars + IGParsAdd
                sd2DiscrFac <- sd2IG <- as.vector(rInvGamma(1, IGPars["alpha"], IGPars["beta"] ))
            })            
    sd2DiscrFac <- mean(sd2DiscrFacStreams)
    ##value<< list with components
    list(	##describe<<
            isUpdated=TRUE		##<< boolean, if changed block parameters, always true with Gibbs sampling
            , xC=log(sd2DiscrFac)	    ##<< numeric vector: components of position in parameter space that are being updated
    )	##end<<
}
.tmp.f <- function(){
    #IGPars <- getInvGammaParsFromMeanAndMode( 20, 0.05 )
    #IGPars <- c(alpha=0.001,beta=0.001)
    xGrid <- seq(0.01, 4, length.out=101)
    plot( (tmpd <- sqrt(dinvgamma(xGrid, IGPars["alpha"], IGPars["beta"]))) ~ xGrid)
    plot( head(tmpd) ~ head(xGrid))
    tmpMean <- IGPars["beta"]/(IGPars["alpha"]-1)
    tmp <- rInvGamma(100, IGPars["alpha"], IGPars["beta"] )
    abline(v=tmpMean, col="grey")
    c(mean(tmp))
    # plotting rInvGamma distribution
}

.depr <- function(){
    #inside .updateLogSigmaDiscrByGibbsSamplingGamma tried to base delta only on data instead of former sigma^2
    # did not work
    #sd2DiscrData <- .estimateSignalVariance( int$zO, int$LambdaOO, int$sd2ObsO)
    sd2DiscrData <- .estimateSignalVarianceOpt( int$zO, int$LambdaOO, int$sd2ObsO)
    KOO <- sd2DiscrData*int$LambdaOO 
    Ky <- KOO + diag(int$sd2ObsO, nrow=length(iO))
    solvedKyZ <- solve(Ky, int$zO)
    deltaOD <- as.vector(KOO %*% solvedKyZ)
    solvedLambdaDeltaO <- sd2DiscrData * solvedKyZ  # only if detlaOD are expected values
    normDeltaO <- deltaOD %*% solvedLambdaDeltaO
}


#' @export
.computeLogDenLogPsi <- function( 
        logPsi, sd2Discr, predProcX, deltaAIntermediate, obs, sd2Obs, priorGammaPsiPars
        #,T=2
        ,weightPSd2Discr  ##<< the weight for signal variance in the likelihood of correlation length psi
){ #, logSd2DiscrVar ){ #, priorPsiPars ){
    logPriorPsi <- as.vector((priorGammaPsiPars["k"]-1)*logPsi -1/priorGammaPsiPars["theta"]*exp(logPsi))
    if( length(attr(deltaAIntermediate,"problems")) ) return(c( 
                        mD = -Inf   
                        ,obs= -Inf
                        ,priorPsi = logPriorPsi
                ))
    deltaA <- deltaAIntermediate$deltaA
    iO <- deltaAIntermediate$iO
    deltaO <- deltaA[ deltaAIntermediate$iO ]
    residuals <- obs - (predProcX + deltaA) 
    #
    logLikObs <- -1/2 * sum( residuals^2/sd2Obs )
    logLikGP <- deltaAIntermediate$logPDeltaO
    res <- c( 
            mD = logLikGP  
            ,obs= logLikObs
            ,priorPsi = logPriorPsi
    )
}


.tmp.f <- function(){
    library(parallel)
    if(exists("cl")) stopCluster(cl)
    cl <- makeCluster(2L)
    #clusterExport(cl, c("mcCommon", ".sampleRangeOneChain"), envir=environment())
    clusterEvalQ(cl, library(blockDemac) ) 
    clusterEvalQ(cl, library(twMisc) ) 
    
    .tmp.f <- function(){
        clusterEvalQ(cl, helloFromBlockDemac() )
        dumpedStop <- dumpOnError(stop)
        lapply( 1, dumpedStop("dumpedStop"))
        load("last.dump.rda")
        debugger()
    }
    clusterExport(cl, c("simpleModel"))
    sampler@cluster <- cl
    sampler <- setupAndSample(sampler, nSample=30L)
    #stopCluster(cl)
    sampleLogs <- getSampleLogs(sampler)
    plot( asMcmc(sampleLogs), smooth=FALSE )
    
}

#' @export
compileLogSd2DiscrBlockDescription <- function(
    ### create Block Descriptions for updating logSd2Discr, the common model discrepancy factor
    discrepancyBlocks   ##<< list of Discrepancy block compilations by \code{\link{compileDiscrepancyBlocksForStream}}        
    ,priorInvGammaSd2DiscrPars = getInvGammaParsFromMeanAndMode( 20, 0.05 ) ##<< prior for corellation length
){
    ##details<< 
    ## the entries of discrepancyBlocks must each contain an element \code{meanSdObs}
    ## giving the sqare root of mean observation variance. Its name must correspond to 
    ## the name of the intermediate providing the model discrepancies deltaO.
    # Gibbs IG block for sd2Discr
    sd2ObsClosure <- unlist(structure(lapply( discrepancyBlocks, "[[", "meanSd2Obs" ),names=NULL))
    namesIntermediatesClosure <- names(sd2ObsClosure)
    logSd2DiscrVar <- "logSd2Discr"
    logSd2DiscrSpec <- blockSpec(logSd2DiscrVar,logSd2DiscrVar,   
            new("FunctionBasedBlockUpdater"
                    , fUpdateBlock=function(...){.updateLogSigmaDiscrByGibbsSamplingGamma(...)}  # allows tracing called function   
                    , argsFUpdateBlock=list( priorInvGammaPars=priorInvGammaSd2DiscrPars
                        , intermediateName=namesIntermediatesClosure
                        , sd2Obs=sd2ObsClosure )  
            )        
            ,intermediatesUsed=namesIntermediatesClosure
    )
}

computeInverseOfBlockedMatrix <- function(
    ### invert a blocked matrix [A,B;C,D]
    A,B,C,D, invA=solve(A)
){
    # see: https://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion
    invE <- D - C %*% invA %*% B
    E <- solve(invE)
    BE <- B %*% E
    CinvA <- C %*% invA
    ##value<< the blocks of the inverted matrix
    list(
            A = invA + invA %*% BE %*% CinvA
            ,B = -invA %*% BE
            ,C = -E %*% CinvA
            ,D = E
        )
}
attr(computeInverseOfBlockedMatrix,"ex") <- function(){
    A <- matrix(c(1,.1,.1,1), nrow=2L)
    B <- matrix(0.05, nrow=2L, ncol=3L)
    C <- t(B)
    D <- matrix(0.1, nrow=3L, ncol=3L); diag(D) <- 1L
    M <- rbind( cbind(A,B), cbind(C,D) )
    
    invM1 <- solve(M)
    invM1 %*% M - diag(nrow=5)
    res <- invBlockedMatrix(A,B,C,D)
    invM2 <- rbind( cbind(res$A,res$B), cbind(res$C,res$D) )
    invM2 %*% M - diag(nrow=5)  # all near 0
    
    b1 <- 1:2
    b2 <- 3:5
    b <- c(b1,b2)
    b %*% M %*% b
    blockedMatrixNorm(b1,b2, list(A=A,B=B,C=C,D=D)) - as.vector(b %*% M %*% b)  # zero difference
}

blockedMatrixNorm <- function(
        ### compute (b1, b2T)^T M (b1, b2) for M of form [A,B;C,D]
        b1,  b2
        , blocks    ##<< list with four entries corresponding to blocks of matrix (by row) [A,B;C,D]
){
    b1T <- t(b1)
    b2T <- t(b2)
    as.vector((b1T %*% blocks[[1]] %*% b1) + (b2T %*% blocks[[3]] %*% b1) + (b1T %*% blocks[[2]] %*% b2)  + (b2T %*% blocks[[4]] %*% b2))
}





