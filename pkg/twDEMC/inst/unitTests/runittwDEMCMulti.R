t_est.ofMultiIntermediateCor <- function(){
    int.ofMultiIntermediate( loadBalancing=FALSE, useConditionalProposal=TRUE )
}


int.ofMultiCor <- function(loadBalancing=FALSE, useConditionalProposal=FALSE){
    require(mvtnorm)
    sigma <-  matrix(c(8,7.9,7.9,8), ncol=2)
    sigma0 <- matrix( c(9,4,4,9), ncol=2 )      # initially broader distributed but already correlated
    thetaMean <- c(a=0,b=0)
    x <- rmvnorm(n=1000, mean=thetaMean, sigma=sigma)
    #plot( x )
    #points(x, col="blue")
    #
    denMVCor <- function(popt, sigma){
        dmvnorm( popt, sigma=sigma, log=TRUE )
    }
    #hist(apply( x, 1, denMVCor, sigma=sigma ))
    argsFLogDen <- list(
            sigma=sigma
    )            
    #
    .nPop=2
    .nChainPop=4
    ZinitPops <- initZtwDEMCNormal( c(a=0,b=0), sigma0, nChainPop=.nChainPop, nPop=.nPop, m0FiniteFac = 0.2)
    #
    .nGen=1024
    #undebug(.sampleStatesConditional)
    .remoteDumpfileBasename="tmp/dump_twDEMC_ofMultiCor"
    resa <- concatPops( resBlock <- twDEMCBlock( 
                    ZinitPops
                    , nGen=.nGen 
                    ,dInfos=list(
                            d1=list(fLogDen=denMVCor, argsFLogDen=argsFLogDen)
                    )
                    ,blocks = list(
                            a=list(dInfoPos="d1", compPos="a")
                            ,b=list(dInfoPos="d1", compPos="b")
                    )
                    ,nPop=.nPop
                    #,controlTwDEMC=list(thin=8, useConditionalProposal=FALSE)		
                    ,controlTwDEMC=list(thin=8, useConditionalProposal=TRUE)		
                    #,debugSequential=TRUE
                    ,remoteDumpfileBasename = .remoteDumpfileBasename
                    ,doRecordProposals=TRUE
            ))
    plot( as.mcmc.list(resa), smooth=FALSE)
    rescoda <- as.mcmc.list(resaE <- thin(resa,start=.nGen%/%2))
    (tmp <- cov(stackChains(resaE)[,-1]))
    #plot(stackChains(resaE)[,-1])
    gelman.diag(rescoda)
    #plot(rescoda, smooth=FALSE)
    #
    #i=1
    (.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
    (.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
    # 1/2 orders of magnitude around prescribed sd for theta
    # no true sdTheta
    #.pop=1
    for( .pop in seq(along.with=.popsd) ){
        checkMagnitude( sqrt(diag(sigma)), .popsd[[.pop]] )
    }
    #
    # check that thetaTrue is in 95% interval 
    .pthetaTrue <- sapply(1:2, function(.pop){
                pnorm(thetaMean, mean=.popmean[[.pop]], sd=.popsd[[.pop]])
            })
    checkInterval( .pthetaTrue ) 
}

