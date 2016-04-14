.calcMaxConstrainedLogDensity <- function(
    ## calculating the logDensity and adjusting results
    fLogDensity         ##<< logDensity function
    , parameters            ##<< parameter vector (first argument to fLogDen)
    , intermediates         ##<< list of intermediate results provided to the logDensity function
    , logDensityComponents  ##<< numeric vector( nResCompDensity): only NA values need to be recomputed
    , argsFLogDensity   ##<< further named arguments to fLogDen
    , maxLogDensity     ##<< numeric vector (nResComp): maximum logDensity (overfitting control, usually -1/2 nObs)
){
    LpOverfit <- do.call( fLogDensity, c(list(parameters, intermediates, logDensityComponents), argsFLogDensity) )	# evaluate logDen
    #trace(logDenThetaSa, recover)      #untrace(logDenThetaSa)
    #do.call( logDenThetaBlocks, c(list(parameters, intermediates, logDensityComponents), argsFLogDensity) )	# evaluate logDen
    lengthDen <- length(maxLogDensity)
    assert_that( lengthDen==0 || lengthDen==1 || lengthDen==length(LpOverfit) )
    Lp <- if( length(maxLogDensity) ) pmin( LpOverfit, maxLogDensity ) else LpOverfit
}


.metropolisDecision2 <- function(
    logDenCompCurrent    ##<< numeric vector (nResComp) result components for current density
    ,logDenCompProposal  ##<< numeric vector (nResComp) result components for density of proposed state
    ,rExtra              ##<< numeric scalar: log of Snooker updates Metropolis ratio correction factor 
    ,tempResCompC        ##<< numeric vector (nResComp) Temperature that the current metropolis ratio is based on
    ,logDenPBeforeDR=numeric(0)  ##<< numeric vector (nResComp): provide logDensity of original step to apply correction for delayed rejection
    ,iInternalLogDensityComponents=integer(0)##<< indices of logDensityComponents handled internally
){
    #if( any(!is.finite(logDenCompCurrent))) recover() 
    if (!all(is.finite(logDenCompProposal))) return(FALSE)
    ##details<< \describe{\item{Internal Metropolis descision}{
    ## TODO
    ## Argument code{iInternalLogDensityComponents}  denotes the indices among the results of fLogDen that are handled internally  
    ## }}
    posTExt <- setdiff( seq_along(tempResCompC), iInternalLogDensityComponents )		#externally handled components
    nExt <- length(posTExt)
    #Metropolis step with temperature and correction of Snooker update (rExtra)
    logrDS20 <- (logDenCompProposal[posTExt] - logDenCompCurrent[posTExt]) / tempResCompC[posTExt]
    logAlpha2 <- logAlpha20 <- logDenSum <- rExtra + sum(logrDS20)
    isDelayedRejection = length(logDenPBeforeDR) != 0
    if( isDelayedRejection ){
        # TODO check Haario06 if correct
        # correct with first stage DR factor (1-alpha21)/(1-alpha10) with meaning 0:accepted 1:first proposal 2:second proposal
        logAlpha10 <- sum( (logDenPBeforeDR[posTExt]-logDenCompCurrent[posTExt])/tempResCompC[posTExt] ) 
        logrDS21 <- (logDenPBeforeDR[posTExt] - logDenCompProposal[posTExt]) / tempResCompC[posTExt]
        logAlpha21 <- sum(logrDS21)
        logAlpha2 <- suppressWarnings( logAlpha20  +log(1-exp(logAlpha21)) -log(1-exp(logAlpha10)) )	# log and exp may produce NaNs
    } 
    runif1 <- runif(1)  # store in local variable instead of using it directly in expression for better debugging
    isAccepted <-  is.finite(logAlpha2) && (logAlpha2 > log(runif1))
    return(isAccepted)
}	

popApplyTwDEMC <- function( 
        ### Applying a function across all chains of one population for each case.
        x				##<< a matrix with columns chains or array with last dimension chain
        ,nPop			##<< number of populations
        ,FUN			##<< function to apply to population submatrix
        ,...			##<< further arguemtns to FUN
){
    #popApplyTwDEMC
    ##seealso<<   
    ## \code{\link{popMeansTwDEMC}}
    ## \code{\link{subChains.twDEMC}}
    ldim <-length(dim(x)) 
    nChainsP <- dim(x)[ldim] %/% nPop
    resl <- 
            if( ldim==2 ) 
                lapply( (0:(nPop-1)), function(iPop0){
                            iChains <-  iPop0*nChainsP + (1:nChainsP)
                            FUN(x[,iChains ,drop=FALSE],...) 		
                        })
            else if( ldim==3 )
                #iPop0=0
                lapply( (0:(nPop-1)), function(iPop0){
                            iChains <-  iPop0*nChainsP + (1:nChainsP)
                            FUN(x[,,iChains ,drop=FALSE],...) 		
                        })
            else stop("dimension of x must be 2 or 3")
    abind(resl, rev.along=0)
    #resl
    ### array with last dimenstion correponding to population
}


