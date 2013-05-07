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
                            dSparse=list(fLogDen=denSparsePrior, argsFLogDen=argsFLogDen, intermediate="modTwoDenEx")
                            ,dRich=list(fLogDen=denRichPrior, argsFLogDen=argsFLogDen, intermediate="modTwoDenEx")
                    )
                    ,blocks = list(
                            a=list(dInfoPos="dSparse", compPos="a")
                            ,b=list(dInfoPos="dRich", compPos="b")
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
