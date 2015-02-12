test.simple <- function(){
    data(twLinreg1)
    
# collect all the arguments to the logDensity in a list (except the first argument of changing parameters)    
    argsFLogDen <- with( twLinreg1, list(
                    fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
                    obs=obs,			    ### vector of data to compare with
                    invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
                    thetaPrior = thetaTrue,	### the prior estimate of the parameters
                    invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
                    xval=xval
            ))
    do.call( logDenGaussian, c(list(theta=twLinreg1$theta0),argsFLogDen))
    do.call( logDenGaussian, c(list(theta=twLinreg1$thetaTrue),argsFLogDen))    # slightly largere misfit than nObs/2=15, underestimated sdObs
    
    .nGen=200
    .nPop=2
    mcp <-  twDEMCSA( 
            theta=twLinreg1$theta0, covarTheta=diag(twLinreg1$sdTheta^2)       # to generate initial population
            , nGen=.nGen
            , dInfos=list(den1=list(fLogDen=logDenGaussian, argsFLogDen=argsFLogDen))
            , nPop=.nPop                                        # number of independent populations
            , controlTwDEMC=list(thin=4)                        # see twDEMCBlockInt for all the tuning options
            , ctrlConvergence=list(maxRelTChangeCrit=0.1)       # ok if T changes less than 10% 
            , ctrlT=list(TFix=c(parms=1))                       # do not use increased temperature for priors
            , nObs=c(obs=length(argsFLogDen$obs))               # number of observations used in temperature calculation
    )
#mcp <- twDEMCSA( mcp, nGen=2000) 
    mcp <- twDEMCSA( mcp, nGen=1400)     # continue run (for real do more generations, not converged yet)
    rescoda <- as.mcmc.list(subsetTail(mcp))
    
    .expSdTheta <- twLinreg1$sdTheta
    .expTheta <- twLinreg1$thetaTrue
    suppressWarnings({	
                summary(rescoda)$statistics[,"Mean"]
                summary(rescoda)$statistics[,"SD"]
                .expTheta
                .expSdTheta
                (.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
                (.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
            })
    
    # 1/2 orders of magnitude around prescribed sd for theta
    .pop=1
    for( .pop in seq(along.with=.popsd) ){
        checkMagnitude( .expSdTheta, .popsd[[.pop]] )
    }
    
    # check that thetaTrue is in 95% interval 
    .pthetaTrue <- sapply(1:2, function(.pop){
                pnorm(.expTheta, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
            })
    checkInterval( .pthetaTrue ) 
    
    #---------- decrease towards temperature 1
    (.T <- getCurrentTemp(mcp))
    .nObs <- c(parms=getNParms(mcp), obs=length(mcp$dInfos[[1]]$argsFLogDen$obs))
    calcBaseTemp( .T, .nObs[names(.T)], TFix=c(parms=1) )       # 15% bias inferred
    mcp0 <-  twDEMCSA( 
            theta=mcp
            , nGen=800
            , ctrlT=list(TFix=c(parms=1), TEndFixed=1)                       # do not use increased temperature for priors
    )
    TBase0 <- calcBaseTemp( getCurrentTemp(mcp0), .nObs[names(.T)], TFix=c(parms=1) ) -1
    checkTrue( TBase0 < 1e-4 )
}

test.ofMultiIntermediate <- function(){
    loadBalancing=FALSE
    # same as in test case twDEMC, but here check that parallel calculation of initial logDen
    data(twTwoDenEx1)
    thresholdCovar = 0.3	# the true value
    #thresholdCovar = 0		# the effective model that glosses over this threshold
    thetaMean <- twTwoDenEx1$thetaTrue
    sdTheta <- (thetaMean*0.3)
    #
    argsFLogDen <- list(
            thresholdCovar=thresholdCovar
            , thetaPrior=thetaMean
            , invCovarTheta=sdTheta^2
            , twTwoDenEx=twTwoDenEx1
    )    
    #"parmsSparse" "obsSparse"   "parmsRich"   "obsRich" 
    nObs <- with( argsFLogDen$twTwoDenEx$obsTrue, c(
                    parmsSparse=1,  obsSparse=length(y1)
                    , parmsRich=1, obsRich=length(y2) 
    ))
    #
    .nGen=128
    .nPop=2
    #undebug(denRichPrior)
    .remoteDumpfileBasename=NULL
    #.remoteDumpfileBasename="dump_twDEMCSA_ofMultiIntermediate"
    resa <- concatPops( resBlock <- twDEMCSA( 
                    thetaPrior = thetaMean*1.5
                    , covarTheta= diag(sdTheta^2)
                    , nObs=nObs
                    , nGen=2*.nGen 
                    ,dInfos=list(
                            dSparse=list(fLogDen=denSparsePrior, argsFLogDen=argsFLogDen)
                            ,dRich=list(fLogDen=denRichPrior, argsFLogDen=argsFLogDen)
                    )
                    ,blocks = list(
                            a=list( compPos="a", dInfoPos="dSparse", intermediateId="modTwoDenEx")
                            ,b=list(compPos="b", dInfoPos="dRich",  intermediateId="modTwoDenEx")
                    )
                    ,nPop=.nPop
                    ,controlTwDEMC=list(thin=2, loadBalancing=loadBalancing)
                    ,ctrlBatch=list(nGen0=24)
                    ,ctrlT=list(TFix=c(parmsSparse=1,parmsRich=1))
                    #,debugSequential=TRUE
                    ,remoteDumpfileBasename = .remoteDumpfileBasename
                    ,doRecordProposals=TRUE
            ))
    #resBlock$diagnostics
    resa <- concatPops( resBlock <- twDEMCSA( resBlock, nGen=.nGen))    
    #plot( as.mcmc.list(resa), smooth=FALSE)
    #matplot( resa$temp, type="l" )
    resaEnd <- thin(resa,start=getNGen(resa)%/%2)
    #plot( as.mcmc.list(resaEnd), smooth=FALSE)
    #matplot( resaEnd$temp, type="l" )
    rescoda <- as.mcmc.list(resaEnd)
    #
    #i=1
    (.popmean <- suppressWarnings(lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]})))
    (.popsd <- suppressWarnings(lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]})))
    # 1/2 orders of magnitude around prescribed sd for theta
    # no true sdTheta
    #.pop=1
    #for( .pop in seq(along.with=.popsd) ){
    #    checkMagnitude( sdTheta, .popsd[[.pop]] )
    #}
    #
    # check that thetaTrue is in 95% interval 
    .pthetaTrue <- sapply(1:2, function(.pop){
                pnorm(thetaMean, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
            })
    checkInterval( .pthetaTrue ) 
}
