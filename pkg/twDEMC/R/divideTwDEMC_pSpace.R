calcPSubSpacesVol <- function(
     ### estimate of the integrated density of the subspace by average density in volume
     ssLogDenL      ##<< list (nPop) of vector: the (temperated) logDensity of each record in a 
    , qPop          ##<< the previous quantile of samples within space
    , iPopsSpace    ##<< list (nSpace): of integer vectors: pops in space
    , nSamplePop=sapply(ssLogDenL, length)    ##<< the number of samples in subspace
){
    ##description<<
    ## subspace and population are used interchangeably here.
    #
    #iSpace <- names(iPopsSpace)[1]
    maxLogDen <- sapply( ssLogDenL, max, na.rm=TRUE)
    nSpace <- length(iPopsSpace)
    logMeanDensSubsSpaces <- lapply( names(iPopsSpace), function(iSpace){
        iPops  <- iPopsSpace[[iSpace]]
        #iPop <-iPops[1] 
        lw <- sapply( iPops, function(iPop){ # average the unnormalized densities
            twLogSumExp(ssLogDenL[[iPop]], shiftUpperBound=TRUE)-log(nSamplePop[iPop])         
        })
        lw - max(lw, na.rm=TRUE)			#divide by a constant on exp scale to avoid numerical errors
    })
    logMeanDensSubs <- do.call(c, logMeanDensSubsSpaces )
    #barplot(logMeanDensSubs)
    #
    # estimate proportion of subspaces in the limiting distribution
    # it will be used to sample from the subspaces 
    # iiSpace=iSpaces[1]
    pPops <- numeric( length(qPop))
    #iiSpace <- 1
    for( iiSpace in 1:nSpace){
        # two components: volume quantile and density importance
        p2u <-  exp(logMeanDensSubsSpaces[[iiSpace]])*qPop[ iPopsSpace[[iiSpace]] ] 		
        # normalize to 1 within population
        pPops[ iPopsSpace[[iiSpace]] ] <- p2 <- pmax( .Machine$double.eps, p2u/sum(p2u)	)	 
    }
    ### numeric vector (nPop): the proportion of the subspace within its space
    pPops
}

calcPSubSpacesHME <- function(
        ### estimate of the integrated density of the subspace by harmonic mean estimate (HME)
        ssLogDenL      ##<< list (nPop) of vector: the (temperated) logDensity of each record in a 
        , qPop          ##<< the previous quantile of samples within space
        , iPopsSpace    ##<< list (nSpace): of integer vectors: pops in space
        , nSamplePop=sapply(ssLogDenL, length)    ##<< the number of samples in subspace
        , qOmitLow=0.05 ##<< quantile of low logDensities to omit in integral to stabilize HME estimate
){
    ##description<<
    ## subspace and population are used interchangeably here.
    #
    #iSpace <- names(iPopsSpace)[1]
    nSpace <- length(iPopsSpace)
    pPops <- numeric(length(qPop))
    for( iSpace in  names(iPopsSpace) ){
                iPops  <- iPopsSpace[[iSpace]]
                #iPop <-iPops[1] 
                lhme <- sapply( iPops, function(iPop){ # average the unnormalized densities
                            ssLogDen <- ssLogDenL[[iPop]]
                            ssLogDenUpper <- ssLogDen[ ssLogDen > quantile(ssLogDen, qOmitLow)]
                            log(length(ssLogDenUpper)) -  twLogSumExp(-ssLogDenUpper, shiftLowerBound=TRUE)    # n / sum e^-L     
                        })
                lw <- lhme - twLogSumExp(lhme)        # hme / sum(hme)
                pPops[iPops] <- exp(lw)
            }
    ### numeric vector weights for each population (sum to one within space)
    pPops
}


