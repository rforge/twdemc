# ' @importFrom MCMCpack rinvgamma # ported to this package because of dependency problems 
#' @importFrom mvtnorm rmvnorm


#' @export
compileDiscrepancyBlocksForStream <- function(
        ### compile the intermediates and blocks for GP model error of one stream
        x           ##<< locations of observations
        ,obs        ##<< observations
        ,sd2Obs     ##<< variance of observation uncertainty
        , predProcSpec   ##<< named list of one intermediate specification for model prediction
        , streamSuffix=""	##<< suffix applied to parameters, and specifications to distinguish data streams
        ,priorGammaPsiPars = getGammaParsFromMeanAndVariance(1,0.2) ##<< prior for variance of model discrepancy as a factor of average observation uncertainty
        ,priorInvGammaSd2DiscrPars = getInvGammaParsFromMeanAndMode( 20, 0.05 ) ##<< prior for corellation length
        #,priorInvGammaSd2DiscrPars = c( alpha=0.001, beta=0.001 ) ##<< prior for signal variance
        ,minSd2DiscrFraction = 0 #1/20   ##<< minimum of uncertainty of GP training data as a fraction of model error variance, deprecated: was used to prevent inflating signal variance by  
        ,xPred=numeric(0)  ##<< locations at which to sample model discrepancies, mapped to paramter vector md_1 to md_<nPred>
        ,isSigmaBlockReturned=TRUE  ##<< set to FALSE to omit block for updating sigma_streamSuffix
        ,meanSd2Obs=mean(sd2Obs)    ##<< magnitude of the variance of observation uncertainty
){
    if( !is.list(predProcSpec) || (length(names(predProcSpec)) != 1L) || !is(predProcSpec[[1]],"IntermediateSpecification")) stopDemac(
                "predProcSpec must be a named list with entry of class IntermediateSpecification.")
    if( length(x) != length(obs) ) stopDemac("locations x and observations obs must be of the same length.")
    GPLimits <- .computeGPLimits(x)
    if( length(sd2Obs)==1L ) sd2Obs <- rep(sd2Obs, length(obs) )
    if( length(sd2Obs) != length(obs) ) stopDemac("sd2Obs and observations obs must be of the same length.")
    # names of used intermediates and parameters including suffix
    predProcName <- names(predProcSpec)[1]
    iOName <- paste("iO",streamSuffix,sep="_")
    #iKyName <- paste("iKy",streamSuffix,sep="_")
    iDiscrepancyName <- paste("iPz",streamSuffix,sep="_")
    #deltaAName <- paste("deltaA",streamSuffix,sep="_")
    logPsiVar <- paste("logPsi",streamSuffix,sep="_")
    dsVar <- paste("ds",streamSuffix,sep="_")
    logSigma2Var <- 
            logSd2DiscrVar <- if( isSigmaBlockReturned ){
                paste("logSd2Discr",streamSuffix,sep="_")
            } else {
                "logSd2Discr"
            }
    #
    meanSd2ObsClosure <- meanSd2Obs #mean(sd2Obs)
    xPredVal <- force(xPred)    # force evaluation to assign to closure
    xVal <- force(x)
    iOSpec <- intermediateSpec(
            function(theta, intermediates ){
                .computeGPTrainingLocations(ds=theta[dsVar], psi=exp(theta[logPsiVar]), x=xVal, xPred=xPredVal )
            }
            ,argsUpdateFunction=list()
            ,parameters = c(dsVar, logPsiVar)
    )
    #
    force(priorInvGammaSd2DiscrPars) # will be needed in closure
    iDiscrepancy <- intermediateSpec(
            function(theta, intermediates ){
                pred=intermediates[[ predProcName]][1:length(obs)] # omit process prediction locations
                .computeSolvedKyZ(obs=obs,  pred=pred
                    ,iORes=intermediates[[ iOName]]
                    , sd2Obs=sd2Obs
                    ,meanSd2Obs=meanSd2ObsClosure
                    ,priorInvGammaSd2DiscrPars = priorInvGammaSd2DiscrPars
                    ,x=x
                )
            }
            ,argsUpdateFunction=list()
            ,parameters = c(dsVar, logPsiVar)  # for tracing errors (else only depending on intermediates)
            ,intermediates = c(predProcName, iOName ) 
    )
    #
    intermediateSpecs = c( predProcSpec, list( iOSpec, iDiscrepancy))
    names(intermediateSpecs) <- c(predProcName, iOName, iDiscrepancyName )
    #
    # update block for spacing of supporting locations dS
    dsSpec <- blockSpec(dsVar,c(dsVar, logPsiVar),   
            new("FunctionBasedBlockUpdater"
                    , fUpdateBlock=function(theta, iParms
                            ,lowerParBounds ##<< named numeric vector, lower limits of the parameters to update  
                            ,upperParBounds ##<< named numeric vector, upper limits of the parameters to update  
                            ,intermediates=list() ##<< intermediate result for current state xC, see end of vignette on using two Densities
                        ){
                            .updateDs(
                                dsPrev=theta[dsVar], psi=exp(theta[logPsiVar])
                                ,minSpacing = max(lowerParBounds, GPLimits$minXSpacing)
                                ,maxSpacing = min(upperParBounds, GPLimits$maxXSpacing)
                                )}     
                    , argsFUpdateBlock=list() 
                    )
    )
    #
    # Metropolis block for logPsi
    namesResLogDenPsi <- paste(c("mD","obs","priorPsi"),streamSuffix,sep="_")
    nX <- length(x) 
    force(priorGammaPsiPars) # will be needed in closure
    logDenLogPsiSpec <- blockSpec(logPsiVar,,          
            new("MetropolisBlockUpdater",
                    fLogDensity=function(theta, intermediates, logDensityComponents, obsVal, sd2ObsVal){
                        res <- .computeLogDenLogPsi( logPsi=theta[logPsiVar]
                                , iPzRes=intermediates[[iDiscrepancyName]]
                                , priorGammaPsiPars=priorGammaPsiPars
                        )
                        names(res) <- namesResLogDenPsi # to distinguish from other streams
                        res
                    }
                    ,argsFLogDensity=list(obsVal=obs, sd2ObsVal=sd2Obs)   #priorPsiPars by closure
                    ,logDensityComponentNames = namesResLogDenPsi
            )
            ,intermediatesUsed=c(iDiscrepancyName)   # adds dependencies to parameters used by intermediates
    )
    #
    blockSpecs <- list(dS=dsSpec, logDenLogPsi=logDenLogPsiSpec)
    if( !isTRUE(isSigmaBlockReturned) ) blockSpecs$logSd2Discr <- NULL
    names(blockSpecs) <- paste(names(blockSpecs),streamSuffix,sep="_")
    #
    xiVars <- character(0)
    nPred <- length(xPred) 
    if( nPred ){
        xiVars <- paste("xi",1:nPred,streamSuffix,sep="_")
        xiSpec <- list( xi=blockSpec(xiVars,,          
                new("FunctionBasedBlockUpdater"
                        , fUpdateBlock=function(theta, iParms
                                ,lowerParBounds ##<< named numeric vector, lower limits of the parameters to update  
                                ,upperParBounds ##<< named numeric vector, upper limits of the parameters to update  
                                ,intermediates=list() ##<< intermediate result for current state xC, see end of vignette on using two Densities
                        ){
                            iPz=intermediates[[iDiscrepancyName]]
                            iORes <- intermediates[[iOName]]
                            predP <- intermediates[[ predProcName]][-(1:length(obs))] 
                            .sampleProcessValues(
                                    Ky=iPz$Ky, solvedKyZO=iPz$solvedKyZO
                                    ,LambdaPO=iORes$LambdaPO, LambdaPP=iORes$LambdaPP
                                    ,sd2Discr = iPz$sd2Discr
                                    ,predP = predP
                                    ,xiPrev = theta[iParms]
                            )}     
                        , argsFUpdateBlock=list()
                )
                , intermediatesUsed=c(iDiscrepancyName, iOName, predProcName)
        ))
        names(xiSpec) <- paste("xi",streamSuffix,sep="_") 
        blockSpecs <- c(blockSpecs, xiSpec)
    }
    fConstrainingSubSpace <- .getFConstrainingSubSpace(GPLimits$logPsiMin,GPLimits$logPsiMax, logPsiVar )
    ##value<< list with entries
    list(
            intermediateSpecs=intermediateSpecs ##<< list with intermediateSpecifications predProc_suffix, iO_suffix, and deltaA_suffix  
            ,blockSpecs=blockSpecs              ##<< update block specifications logDenLogPsi_suffix, and logSd2Discr_suffix 
            ,GPLimits=GPLimits                  ##<< results of \code{\link{computeGPLimits}}
            ,fConstrainingSubSpace=fConstrainingSubSpace
            ,namesProcessSamples= xiVars        ##<< names of the parameters of process samples 
            ,meanSdObs=structure( sqrt(meanSd2ObsClosure),names=iOName) ##<< named numeric scalar: average observation uncertainty, with name of deltaA intermediate for this stream 
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
        confint(lm(obs ~x-1))
        theta0 = c(theta=as.numeric(thetaHat), ds=1, logPsi=log(0.6))
        covarTheta0 = diag(c(theta0*0.2)^2, nrow=length(theta0))    # spread for initial population
        # ---- define the prediction intermediate 
        xPred <- c(1.5,3.8,6)
        #xPred <- numeric(0)
        nPred <- length(xPred)
        predProcSpec <- intermediateSpec(
                dumpOnError(function(theta, intermediates, xVal ){
                            simpleModel(theta["theta"],xVal)
                        })
                ,argsUpdateFunction=list(
                        xVal=c(x,xPred)
                    )
                ,parameters = c("theta")
                ,intermediates = character(0)  
        )
        streamSuffix <- "obs1"
        theta0N <- c(theta0, if(nPred) structure(rep(0, length(xPred)),names=paste("xi",seq_along(xPred),sep="_")) else numeric(0) )
        names(theta0N)[2:(3+nPred)] <- paste(names(theta0N)[2:(3+nPred)],streamSuffix,sep="_")
        covarTheta0 =  diag(c(c(0.1, 0.1, 0.1)^2, rep(0.2,nPred))^2, nrow=length(theta0N))    # spread for initial population
        #
        res <- compileDiscrepancyBlocksForStream(x, obs, sd2Obs, predProcSpec = list(predSpeed=predProcSpec), streamSuffix, xPred=xPred)
        #
        predProcRes0 <- getUpdateFunction(res$intermediateSpecs$predSpeed)(theta0, list(), xVal=c(x,xPred)) 
        iORes0 <- getUpdateFunction(res$intermediateSpecs$iO_obs1)(theta0N, list()) 
        iPzRes0 <- getUpdateFunction(res$intermediateSpecs$iPz_obs1)(theta0N, list(iO_obs1=iORes0, predSpeed=predProcRes0))
        intermediates0 <- list(predSpeed=predProcRes0, iO_obs1=iORes0, iPz_obs1=iPzRes0 )
        #
        testLogDenLogPsi <- getFLogDensity(getBlockUpdater(res$blockSpec$logDenLogPsi_obs1))(
                theta0N
                , intermediates0, numeric(1)
                , obsVal=obs, sd2ObsVal=sd2Obs
        )
        if( nPred ){
            testUpdateXi <- getFUpdateBlock(getBlockUpdater(res$blockSpec$xi_obs1))(
                    theta0N
                    , numeric(0),numeric(0)
                    , intermediates=intermediates0
            )
        }
        #
        # Metropolis block sampling theta, based on intermediate iPz
        logDenTheta <- function( theta, intermediates, logDensityComponents, obs, sd2Obs ){
            pred <- intermediates$predSpeed
            iPz <- intermediates$iPz_obs1
            #sd2Obs <- exp(theta["logSd2Obs"])
            return( c(obsS=iPz$logPzO, obsR=iPz$logPzR, sd2=length(obs)*iPz$logPsd2) )
            #return( c(obsS=iPz$logPzO+iPz$logPzR, obsR=length(x)*iPz$logPsd2) )
            #c(obs=-1/2*sum( (pred - obs)^2 / sd2Obs ))
        }
        testLogDenTheta <- logDenTheta(theta0, intermediates0, numeric(1), obs=obs, sd2Obs=sd2Obs )
        logDenThetaSpec <- blockSpec("theta",,    
                new("MetropolisBlockUpdater",
                        fLogDensity=function(...){logDenTheta(...)},
                        argsFLogDensity=list( obs=obs, sd2Obs=sd2Obs),
                        logDensityComponentNames = names(testLogDenTheta)
                )
                ,intermediatesUsed=c("predSpeed","iPz_obs1")
        )
        #
        res <- compileDiscrepancyBlocksForStream(x, obs, sd2Obs, predProcSpec = list(predSpeed=predProcSpec), streamSuffix, xPred=xPred )
        intermediateSpecs <- c(res$intermediateSpecs)
        blockSpecs <- c(list(bTheta=logDenThetaSpec), res$blockSpecs)
        #mtrace(logDenTheta, recover)   #untrace(logDenTheta)
        theta0N["theta"] <- 0.68
        sampler <- newPopulationSampler( blockSpecs, theta0N, covarTheta0
                ,intermediateSpecifications=intermediateSpecs
                , subSpace=res$fConstrainingSubSpace(new("SubSpace"))
        )
        
        sampler <- setupAndSample(sampler, nSample=48L)
        sampler <- setupAndSample(sampler, nSample=8L)
        #sampler <- setupAndSample(sampler, nSample=8L)
        #sampler <- setupAndSample(sampler, nSample=200L)
        sampleLogs <- getSampleLogs(sampler)
        plot( asMcmc(sampleLogs), smooth=FALSE )
        plot( asMcmc(subsetTail(sampleLogs,0.6)), smooth=FALSE )
        stacked <- adrop(getParametersForPopulation(stackChainsInPopulation(subsetTail(sampleLogs, 0.6)),1L),3)
        plot(density(stacked["theta",]))
        quantile( stacked["theta",], probs=c(0.025,0.975))
        tmp <- t(getLogDensityComponentsForPopulation( subsetTail(sampleLogs, 0.6), 1L )[,,1L])
        
        thetaMean <- rowMeans(stacked)
        thetaCf <- apply(stacked, 1, quantile, probs=c(0.025,0.975))
        yPredTrue <- simpleModel(thetaTrue,xPred)/(1+xPred/20)
        if( nPred ){
            yPredSampled <- thetaMean[3+(1:length(xPred))]
            cfYPredSampled <- thetaCf[,3+(1:length(xPred)), drop=FALSE]
        } else {
            yPredSampled <- numeric(0) #thetaMean[3+(1:length(xPred))]
            cfYPredSampled <- numeric(0)
        }
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



.computeGPLimits <- function(
        ### computing limits on estimates of correlation length psi and spacing of trainign points for Gaussian Process
        x                   ##<< locations of observations
        , maxKernelDim=60   ##<< maximum of dimensionality of the Correlation matrix, corresponding to maximum number of supporting locations
        , psiMin=0          ##<< possibility to specify a minimum correlation length to avoid large matrices on large data streams
){
    ##details<< 
    ## For large psi, some matrices become indeterminate, hence bound by range of x-values.
    ## For very small psi, model error becomes zero between supporting points. However, 
    ## we want to have a smooth model error between supporting points, else we model the noise.
    ## Hence, bound psi by only a bit lower than minimla allowed average distance between supporting points.
    ## If every observation becomes a supporting point, then we model the noise.
    ## Hence, we ensure not allow smaller distance between supporting points than 4 times the average distance
    ## covered by the observations.
    xRange <- diff(range(x))
    maxXSpacing <- xRange / (3+2) 
    minXSpacing = max( psiMin*1.5, xRange/maxKernelDim, min( maxXSpacing, xRange/length(x)*4))
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
.computeGPTrainingLocations <- function(
        ### compute subsets of locations, xO, where observation error will be sampled from Gaussian Process
        ds             ##<< distance between supporting locations
        ,psi           ##<< numeric scalar: correlation length 
        ,x               ##<< locations, e.g. time, at witch observations are available
        ,xPred=numeric(0)##<< location for which to draw process samples (model + delta)
){
    # if suggested psi that is much larger than spacing, this will result in singularities
    if( psi/ds > 6 ) return(list(
                        iO = NA_real_
                        ,LambdaOO = NA_real_
                        ,LambdaRO = NA_real_
                        ,LambdaPO = NA_real_
    ))     
    # determine subset of xO by taking points all every ~psi locations
    xd <- c(0,cumsum(diff(x)))
    iO <- c(1, 1+which( diff(signif(xd,6) %/% signif(ds,6)) != 0))
    # include the edges, else will tend towards zero
    nX <- length(x)
    mIO <- max(iO)
    if( x[nX] > (x[mIO]+ds*3/4) ) iO <- c(iO,nX) else iO[length(iO)] <- nX
    xO <- x[iO]
    LambdaOO <- outer(xO, xO, .corG, psi=psi) 
    LambdaRO <- outer(x[-iO], xO, .corG, psi=psi)
    if( length(xPred) ){
        LambdaPO <- outer(xPred, xO, .corG, psi=psi)
        LambdaPP <- outer(xPred, xPred, .corG, psi=psi)
    }else LambdaPO <- LambdaPP <- numeric(0)
    ##value<< list with entries     
    list(
            iO = iO                 ##<< indices within x, for which model discrepancies are sampled
            ,LambdaOO = LambdaOO
            ,LambdaRO = LambdaRO
            ,LambdaPO = LambdaPO   ##<< correlations among c(xPred,x)
            ,LambdaPP = LambdaPP   
    )
}

# need to export to work on cluster
#' @export
.corG <- function(x1,x2,psi=psi0){
    exp(-((x1-x2)/psi)^2)
}


# need to export to work on cluster
#' @export
.computeSolvedKyZ <- function(
        ### compute Ky^(-1) z
        obs       ##<< numeric vector: observations 
        ,pred     ##<< numeric vector: prediction at c(x,xPred) 
        ,iORes     ##<< index in obs and pred of supporting locations
        ,sd2Obs      ##<< numeric vector (length(obs)) of observation variance 
        ,meanSd2Obs  ##<< numeric scalar: average observation variance
        ,priorInvGammaSd2DiscrPars  ##<< parameters of the prior of the signal variance
        ,x      ##<< locations of obs and pred (for debugging)
){
    z <- (obs - pred)
    iO <- iORes$iO
    nO <- length(iO)
    NAResult <- list(
            z=z
            ,sd2Discr = NA_real_
            ,Ky=NA_real_
            ,solvedKyZO = NA_real_
            ,deltaO = NA_real_
            ,deltaR = NA_real_
            ,logPzO = -Inf
            ,logPzR = -Inf
            ,logPzOK = -Inf
            ,logPsd2 = -Inf
    )
    # if no spacing was calculated return NA
    if( !is.finite(iO[1]) ) return(NAResult)
    zO <- z[iO]
    LambdaOO <- iORes$LambdaOO
    # if singular, return NA
    sd2Discr <- try(.estimateSignalVariance(zO, LambdaOO), silent=TRUE)
    if( inherits(sd2Discr, "try-error") ) return(NAResult)
    #
    KOO <- sd2Discr * LambdaOO
    Ky <- KOO + diag(sd2Obs[iO], nrow=nrow(LambdaOO))
    #
    solvedKy  <- solve(Ky, cbind(zO, KOO))
    solvedKyZO  <- solvedKy[,1]
    solvedKyKOO <- solvedKy[,1+1:nO]
    #deltaO <- as.vector(KOO %*% solvedKyZO)
    #-- really sample deltaO from distribution
    deltaOMu <- as.vector(KOO %*% solvedKyZO)
    deltaOCov <- KOO - KOO %*% solvedKyKOO
    deltaOCovSymm <- (deltaOCov + t(deltaOCov)) / 2     # else sometimes rmvnorm does complain 
    deltaO <- pred <- as.vector(rmvnorm(1, deltaOMu, deltaOCovSymm))   # realization
    #plot(deltaO ~ deltaOMu)
    zR <- z[-iO]
    if( length(zR)) {        
        KRO <- sd2Discr * iORes$LambdaRO
        #deltaR <- as.vector(KRO %*% solvedKyZO)
        #set expected on sampled deltaO instead of expected value given data
        # here do not account for observation uncertainty (use KOO and instead of Ky), because using deltaO instead of y 
        deltaR <- KRO %*% solve(KOO, deltaO)
    } else {
        deltaR <- numeric(0)
    }
    #plot( c(deltaO,deltaR)[order(c(x[iO],x[-iO]))] ~ x, type="l", ylim=range(deltaO,deltaR,z) ); points(z ~ x, pch="+"); abline(h=0, col="grey");  points(deltaOMu ~ x[iO], col="maroon" );  points(deltaO ~ x[iO], col="green", pch="o");      points(deltaR~x[-iO],pch="o");
    #
    detKy <- det(Ky)
    logPzO <- -1/2*( t(zO) %*% solvedKyZO )
    #logPzO <- -1/2*sum( (zO - deltaO)^2/sd2Obs[iO] )
    logPzOK <- logPzO -1/2*(log(detKy) + length(zO)*log(2*pi))
    # non-supporting positions -> Gaussian of process residuals
    logPzR <- if( !length(zR)) 0.0  else {        
                logPzR <- -1/2*sum( (zR - deltaR)^2/sd2Obs[-iO] )  
            }
    # penalty on signal variance
    expectedSd2DiscrFac <- sd2Discr / meanSd2Obs
    logPriorInvGammaSd2Discr <- as.vector(log(expectedSd2DiscrFac)*(-priorInvGammaSd2DiscrPars["alpha"]-1) - priorInvGammaSd2DiscrPars["beta"]/expectedSd2DiscrFac)
    #logPriorInvGammaSd2Discr <- 0
    #
    list(
            z=z
            ,sd2Discr = sd2Discr
            ,Ky=Ky
            ,solvedKyZO = solvedKyZO
            ,deltaO = deltaO
            ,deltaR = deltaR
            ,logPzO = logPzO         ##<< logDen of z given GP and nonchaning n or K
            ,logPzR = logPzR        ##<< logDen of normal of process residual of remaining locations
            ,logPzOK = logPzOK      ##<< logDen of z given GP including changing K and number of supporting locations
            ,logPsd2 = logPriorInvGammaSd2Discr     ##<< prior of observed signal variance
    )
}

.estimateSignalVariance <- function(
        ### estimate signal variance for given correlations and data
        zO      ##<< numeric vector: model data misfit 
        , LambdaOO=outer(xO, xO, .corG, psi=psi) ##<< correlation matrix  
        , psi   ##<< alternatively to specifying LambdaOO, correlation length psi can be specified together with xO or spacing ds
        , xO=.computeGPTrainingLocations(ds, psi, x)    ##<< supporting locations of the GP  
        , ds    ##<< distance between supporting locations
        , x     ##<< locations of observations
){
    solvedLambdaOOZO <- solve(LambdaOO, zO)
    sd2Discr <- as.vector(zO %*% solvedLambdaOOZO) / length(zO)
}

.computeDiscrepancyVariance <- function(
        ### compute variance of model discrepancy 
        parameters  ##<< numeric matrix with columns theta, psi, and ds
        ,fModel     ##<< function of parameter vector, and locations to get model prediction 
        ,obs        ##<< observations
        ,x          ##<< location of observations
        ,colsTheta   ##<< character or integer vector specifying the columns of model parameters
        ,colDs=paste("ds",streamSuffix,sep="_") ##<< column of the distance of supporting locations
        ,colLogPsi=paste("logPsi",streamSuffix,sep="_") ##<< column of the logarithm of correlation length
        ,streamSuffix="obs1"	##<< character scalar
){
    parsStream <- cbind( parameters[,colsTheta], ds=parameters[,colDs], psi=exp(parameters[,colLogPsi]) )
    thetaAll <- parsStream[1,]
    apply( parsStream, 1, function(thetaAll){
        resIO <- .computeGPTrainingLocations(ds=thetaAll["ds"], psi=thetaAll["psi"], x=x)
        iO <- resIO$iO
        LambdaOO <- resIO$LambdaOO
        xO <- x[iO]
        theta <- thetaAll[colsTheta]
        predO <- fModel(theta)[iO]
        zO <- obs[iO] - predO
        sd2 <- .estimateSignalVariance(zO, LambdaOO, thetaAll["psi"], xO, thetaAll["ds"], x)
        sd2
    })
}

.tmp.f <- function(){
    # assume sample1 defined
    fModSparse <- function(theta){
        #trace(modTwTwoDenEx1, recover)  #untrace(modTwTwoDenEx1)
        tmp <- modTwTwoDenEx1(theta, xSparse=xSparse, xRich=xRich)
        tmp$y1
    }
    sampleSd2Sparse <- .computeDiscrepancyVariance(sample1, fModel=fModSparse, obs=obs$y1,x=xSparse, colsTheta=1:2, streamSuffix="sparse" )
    #
    fModRich <- function(theta){
        #trace(modTwTwoDenEx1, recover)  #untrace(modTwTwoDenEx1)
        tmp <- modTwTwoDenEx1(theta, xSparse=xSparse, xRich=xRich)
        tmp$y2
    }
    sampleSd2Rich <- .computeDiscrepancyVariance(sample1, fModel=fModRich, obs=obs$y2,x=xRich, colsTheta=1:2, streamSuffix="rich" )
    plot( density(sampleSd2Sparse/sdObs$y1) )
    lines( density(sampleSd2Rich/sdObs$y2), col="blue" )
    plot( density(log(sampleSd2Sparse/sdObs$y1)) )
    lines( density(log(sampleSd2Rich/sdObs$y2)), col="blue" )
}            




.tmp.f <- function(){
    t <- seq_along(z)
    plot( z ~ t)
    points( iSolvedKyZRes$deltaR ~ t[-iO], col="blue" )
    points( iSolvedKyZRes$deltaO ~ t[iO], col="red" )
    
    plot(obs ~ t)
    points( obs[iO] ~ t[iO], col="red")
    points( pred[iO]+iSolvedKyZRes$deltaO ~ t[iO], col="orange")
    points( pred[-iO]+iSolvedKyZRes$deltaR ~ t[-iO], col="lightblue")
}

#' @export
.updateDs <- function( 
        ### update spacing of supporting locations dS
        dsPrev              ##<< previous value of spacing 
        ,psi                ##<< correlation length to which spacing should adapt
        ,minSpacing         ##<< minimum spacing length
        ,maxSpacing         ##<< maximum spacing length
        , psiSpacingFac=1.5  ##<< multiple of psi, of xSpacing
        , relaxationFac=0.2  ##<< magnitude of jump to target (smaller: slower adaptation) 
){
    ##details<< Target of spacing is 1.5 times psi.
    ## However, let the spacing only slowly apdat to more rapidely fluctuating psi.
    ## If only change of less than 10% is required than do not update spacing
#debug: if( (psi > exp(-6)) && (minSpacing < 0.16) ) recover()    
    dsTarget <- psiSpacingFac * psi
    dsTarget <- min(maxSpacing, max(minSpacing, dsTarget))  # limit to allowed range
    dsChange <- relaxationFac*(dsTarget-dsPrev)
    #XX TODO remove no change when finished
    #dsChange <- 0.0
    relChange <- dsChange/dsPrev
    ##value<< list with components
    #if( relChange > -0.05 && relChange < 0.05) list (
    if( FALSE ) list (
                                    # on small changes do not update ds
                isUpdated=FALSE		        
                , xC=dsPrev	    
        ) else list(	
                isUpdated=TRUE		        ##<< boolean, if changed block parameters, always true with Gibbs sampling
                , xC=dsPrev + dsChange	    ##<< numeric vector: components of position in parameter space that are being updated
        )	
}

#' @export
.sampleProcessValues <- function( 
        ### sample process values (model + discrepancy)
        Ky              ##<< correlation matrix
        ,solvedKyZO     ##<< solve(Ky, zO)
        ,LambdaPO       ##<< correlation matrix of prediction locations with supporting locations
        ,LambdaPP       ##<< correlation matrix among prediction locations 
        ,sd2Discr       ##<< signal variance (unnormalized)
        ,predP          ##<< model predictions for process prediction locations
        ,xiPrev         ##<< previous process samples (returned if not updated)
){
    if( !is.finite(sd2Discr) ) return(list(
            isUpdated=TRUE
            ,xC=xiPrev
       ))
    KPO <- sd2Discr*LambdaPO
    KPP <- sd2Discr*LambdaPP
    muS <- as.vector(KPO %*% solvedKyZO)
    sigmaS <- KPP - KPO %*% solve(Ky, t(KPO))
    xi <- as.vector(rmvnorm(1L, muS, sigmaS))
    list(
            isUpdated=TRUE
            ,xC=predP + xi
            )
}

#' @export
.computeLogDenLogPsi <- function(
        logPsi              ##<< value of psi to calculate logLikelihood for
        ,priorGammaPsiPars  ##<< numeric vector of parameters of prior distribution      
        ,iPzRes     ##<< result of .computePz with entries logPzO, and logPzR

){
     res <- c( 
            mD = iPzRes$logPzOK 
            ,obs= iPzRes$logPzR
            ,priorPsi = as.vector((priorGammaPsiPars["k"]-1)*logPsi -1/priorGammaPsiPars["theta"]*exp(logPsi)) 
    )
}







