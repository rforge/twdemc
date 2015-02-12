orderLogDen <- function(
    ### rank the log-Densities, starting with highest densities
    ss          ##<< numeric matrix ( nRec x nDen+nPar ) 
    ,nDen=attr(ss,"nBlock")     ##<< number of densities
    ,decreasing=TRUE            ##<< argument to order, to start with the highest ranks, i.e best models
){
    ##seealso<<
    ## \code{\link{getBestModelIndices}}
    ## \code{\link{sumLogDenCompBlocks}}
    #        
    oR <- if( nDen==1){
        order( ss[,1], decreasing=decreasing )  
    }else{
        r <- apply( ss[,1:nDen], 2, rank )
        minR <- apply(r,1,min)
        order(minR, decreasing=decreasing)
    }
    ### numeric vector of indices of records with lowest rank of log-Density starting with the maximum density.
    ### For each record the maximum of the density ranks is calculated
    oR
}
attr(orderLogDen,"ex") <- function(){
    if( FALSE ){    # no twDEMCSA result available for now
        #assume resBlock is a result of twDEMCSA with multiple densities
        ss <- stackChains(resBlock, useTemperatedLogDen=TRUE)         # get the temperated summed densities for each block
        ss[ head(orderLogDen(ss)), ]     # display the best results
    }
}


getBestModelIndices <- function(
        ### select the best models based on (temperated) logDensity components
        resLogDen	##<< numeric matrix (nStep x nResComp): logDensity (highest are best)
        , dInfos 	##<< list of lists with entry resCompPos (integer vector) specifying the position of result components for each density
        , prob=0.1  ##<< proportion of good indices to return
        , doOrder=FALSE    ##<< if set to TRUE results will be ordered, yield a slight performance penalty
){
    ##seealso<<
    ## \code{\link{orderLogDen}}
    iBest <- if( length(dInfos) > 1){
                # with several densities, each parameter vector is ranked differently
                # select the case where the maximum of all the ranks across densities is lowest
                #iDen = 1
                logDenBlocks <- sumLogDenCompBlocks(resLogDen, dInfos=dInfos)
                rankDen <- apply( -logDenBlocks,2, rank) # starting with the highest density
                rankDenMax <- apply( rankDen, 1, max )   
                iBest <- which( rankDenMax <= quantile( rankDenMax, prob))
                iBestO <- if( doOrder ) iBest[ order(rankDenMax[iBest]) ] else iBest
            }else{
                logDen <- rowSums(resLogDen)
                iBest <- which( logDen >= quantile( logDen, (1-prob) ))
                iBestO <- if( doOrder ) iBest[ order(logDen[iBest], decreasing=TRUE) ] else iBest
            }
    ### the indices within resLogDen with most highly ranked models (i.e whose worst rank across density is highest)
    iBest
}
attr(getBestModelIndices,"ex") <- function(){
    logDenT <- cbind( -sample(5)/2, -sample(5), -sample(5) )
    #dInfos <- list( d1=list(resCompPos=1:2), d2=list(resCompPos=3) )
    dInfos <- list( d1=list(resCompPos=2), d2=list(resCompPos=3) )
    getBestModelIndices(logDenT, dInfos, 0.2)
}

getBestModelIndex <- function(
        ### select the best model based on (temperated) logDensity components
        logDenT		##<< numeric matrix (nStep x nResComp): logDensity (highest are best)
        , dInfos 	##<< list of lists with entry resCompPos (integer vector) specifying the position of result components for each density 
){
    ##seealso<< 
    ## \code{\link{orderLogDen}}
    iBest <- if( length(dInfos) > 1){
                # with several densities, each parameter vector is ranked differently
                # select the case where the maximum of all the ranks across densities is lowest
                logDenBlocks <- sumLogDenCompBlocks(logDenT, dInfos=dInfos)
                rankDen <- apply( -logDenBlocks,2, rank) # starting with the highest density
                rankDenMax <- apply( rankDen, 1, max )   
                iBest <- which.min(rankDenMax)
            }else{
                iBest <- which.max( rowSums(logDenT) )
            }
    ### the index within logDenT with the best rank
    iBest
}
attr(getBestModelIndex,"ex") <- function(){
    logDenT <- cbind( -sample(5)/2, -sample(5), -sample(5) )
    #dInfos <- list( d1=list(resCompPos=1:2), d2=list(resCompPos=3) )
    dInfos <- list( d1=list(resCompPos=2), d2=list(resCompPos=3) )
    getBestModelIndex(logDenT, dInfos)
    -logDenT
}




sumLogDenCompBlocks <-  function(
        ### sum logDensity components across blocks
        resLogDen	##<< numeric matrix (nStep x nResComp): logDensity (highest are best)
        , dInfos 	##<< list of lists with entry resCompPos (integer vector) specifying the position of result components for each density
){
    ##seealso<<
    ## \code{\link{orderLogDen}}, \code{\link{calcTemperatedLogDen}}
    logDenBlocks <- if( length(dInfos) > 1){
                #iDen = 1
                logDenBlocks <- do.call( cbind, lapply( seq_along(dInfos), function(iDen){
                                    dInfo <- dInfos[[iDen]]
                                    logDenBlock <- rowSums(resLogDen[,dInfo$resCompPos ,drop=FALSE])	
                                }))
            }else{
                logDenBlocks <- matrix(rowSums(resLogDen), ncol=1)
            }
    logDenBlocks
}


.parBoundsEnvelope <- function(
        ### get the parameter bounds that encompasses all given parBounds 
        popsParBounds	##<< list of populations with entries upperParBounds and upperParBounds both named numeric vectors 
){
    #stop(".parBoundsEnvelope: not implemented yet.")
    #parName <- colnames(pop1$parms)[1]
    ub <- lb <- list()
    #pop <- popsParBounds[[2]]
    pnames <- unique( do.call(c, c( 
                            lapply(popsParBounds, function(pop){ names(pop$upperParBounds) }) 
                            ,lapply(popsParBounds, function(pop){ names(pop$lowerParBounds) }) 
                    )))
    #parName <- pnames[1]
    for( parName in pnames ){
        ubPop <- lapply( popsParBounds, function(pop){ pop$upperParBounds[parName] })
        lbPop <- lapply( popsParBounds, function(pop){ pop$lowerParBounds[parName] })
        if( all(sapply(ubPop, function(x){ !is.null(x) && is.finite(x)})) ) 
            ub[parName] <- max(unlist(ubPop)) 
        if( all(sapply(lbPop, function(x){ !is.null(x) && is.finite(x)})) ) 
            lb[parName] <- min(unlist(lbPop)) 
    }
    ##value<< list with entries
    list( ##describe<<
            upperParBounds = unlist(ub)		##<< named numeric vector of upper parameter bounds
            ,lowerParBounds = unlist(lb)	##<< named numeric vector of lower parameter bounds
    ) ##end<<
}
#attr(.parBoundsEnvelope,"ex") <- function(){
.tmp.f <- function(){    
    data(den2dCorEx)
    mc0 <- den2dCorEx$mcSubspaces0
    #mtrace(.parBoundsEnvelope)
    .parBoundsEnvelope( mc0$pops[1:4] )
    .parBoundsEnvelope( mc0$pops[3:4] )
}


