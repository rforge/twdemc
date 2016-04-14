#' @include JumpProposer.R
#' @include DESettings.R
#' @include SampleLogs.R
#' @include AcceptanceTracker.R

setClass("DEJumpProposer", contains="JumpProposer"
        ,representation(  
                deSettings="DESettings"
                ,iParametersThatRequireJumpProposal="integer"
                ,isBurnin="logical"
                ,is1DParameterVector="logical"	    # TRUE if nPar=1, with potential drop of dimensions in slicing
                ,gammaAdjustmentPops="numeric"      # adjustment of jumping length for each population controlled by acceptance rate
                ,acceptanceTracker="AcceptanceTracker"
        )
        ,prototype(isBurnin=TRUE, is1DParameterVector=FALSE, gammaAdjustmentChains=1)
)

setMethod("adjustToAcceptanceRate", signature=c("DEJumpProposer"), 
        function(object,
                ### compile jumping steps
                acceptanceTracker       ##<< object to calculate acceptance rates 
        ) {
            object@acceptanceTracker <- acceptanceTracker
            if( isTRUE(object@deSettings["isAdjustingStepLengthToAcceptanceRate"])) object <- .adjustGammaByAcceptanceRate(object)
            object
        }
)

.adjustGammaByAcceptanceRate <- function(
        ### adjust jumping length Gamma 
        object      ##<< DEJumpProposer object
){
    ##details<<
    ## If acceptance rate is too low, decrease jumping length
    ## If acceptance rate recovers, increase jumping lengh again
    # object@gammaAdjustmentPops will be multiplied to jumping length gamma
    popAcceptanceRates <- computePopulationAcceptanceRates(object@acceptanceTracker)
    iLow <- which( popAcceptanceRates < 0.15 & object@gammaAdjustmentPops > 0.1)
    iHigh <- which( popAcceptanceRates > 0.35 & object@gammaAdjustmentPops < 1.5)
    object@gammaAdjustmentPops[iLow] <- object@gammaAdjustmentPops[iLow] * 0.9 
    object@gammaAdjustmentPops[iHigh] <- object@gammaAdjustmentPops[iHigh] * 1.11 
    object
} 



setMethod("proposeJumps", signature=c("DEJumpProposer"), 
        function(object,
                ### compile jumping steps
                nGeneration
                ,sampleLogs              ##<< SampleLogs object
                ,iPopulationsStep        ##<< integer vector of population indices taking part in this step
                ,iCurrentSamples         ##<< index of the current sample 
        ) {
            if( !length(object@iParametersThatRequireJumpProposal) ){
                sDim <- getSampleDimensions(sampleLogs)
                nChain <- getNChainInPopulation(sDim)*length(iPopulationsStep)
                resNull <- list(jump=array( 0, dim=c(getNParameter(sDim),nGeneration,nChain), dimnames=list(getParameterNames(sDim),NULL,NULL))
                    ,logSnookerDensityMultiplier=matrix(0, nrow=nGeneration, ncol=nChain)
                )
                return(resNull)
            }
            .sampleJumpsFromHistory( object, nGeneration, sampleLogs, iPopulationsStep=iPopulationsStep, iCurrentSamples)
        }
)

setMethod("initializeDimensionalSettings", signature=c("DEJumpProposer"), 
        function(object
                , sampleDimensions                      ##<< class SampleDimension
                , thin=object@deSettings                                  ##<< thinning inteveral, used in calculation of recent history                
                , iParametersThatRequireJumpProposal    ##<< integer vector: indices within parameter vector of parameters for which jumps need to be prepared
        ){
            if( missing(sampleDimensions) ) stop("Missing argument sampleDimensions")
            if( !is(sampleDimensions,"SampleDimensions")) stop("Argument sampleDimensions must be of class SampleDimensions.")
            object@deSettings <- initializeDimensionalSettings(object@deSettings
                    ,nParameterWithProposal = getNParameterWithProposal(sampleDimensions)
                    ,nChainInChainGroup = getNChainInPopulation(sampleDimensions) %/% getNChainGroupsInPopulation(sampleDimensions)
                    ,thin = thin
            )
            object@is1DParameterVector <- (1L == getNParameter(sampleDimensions))
            object@iParametersThatRequireJumpProposal <- iParametersThatRequireJumpProposal
            object@gammaAdjustmentPops <- rep(1L, getNPopulation(sampleDimensions))
            object
        })

if(!exists(".getNStepBack")) setGeneric(".getNStepBack", function(object, acceptanceRates,...) standardGeneric(".getNStepBack"))
setMethod(".getNStepBack", signature=c("DEJumpProposer", "numeric"),
        function(object, acceptanceRates){
            ##details<< 
            ## terBraak report that may not sample the distribution, if not using the full past
            ## but together with decreasing temperature acceptance rate drops very low if 
            ## using the full past
            ## hence constrain to recent past during burnin.
            ## Here, sample among at least nRepresentativeRows recored steps, i.e. nRepresentativeRows*thin generations. 
            ## If Acceptance rate drops then sample among more states up to double nRepresentativeRows
            if( object@isBurnin ){
                ceiling( object@deSettings["nRepresentativeRows"] / pmin(1,acceptanceRates+0.5))
            } else {
                rep( .Machine$integer.max, length(acceptanceRates) )
            }
        })

if(!exists(".sampleJumpsFromHistory")) setGeneric(".sampleJumpsFromHistory", function(object, ...) standardGeneric(".sampleJumpsFromHistory"))
setMethod(".sampleJumpsFromHistory", signature=c("DEJumpProposer"),
        function(object
                ,nGeneration
                ,sampleLogs              ##<< SampleLogs object
                ,iPopulationsStep        ##<< integer vector of population indices taking part in this step
                ,iCurrentSamples         ##<< integer vector (nPopulations): index of the current sample 
        ){
            # -- sample random vectors (note here parameters in rows) 
            # Calculate proposed steps (jumps, i.e. differences instead of absolute destinations) 
            # Done for the entire next thinning interval, instead of only next step
            # Becaus better performance of several step in remote process,
            # and the history of states for generating differences does not change this much
            # across thinningn interval.
            # 
            # zx for ranomd states per step sampled within groups of chains:
            #   random states (nParm,(steps*nChainPop), 4)
            #   first dimension is the state vector
            #   second dimension is steps
            #   third dimension is chains
            ##details<<
            ## make sure to update acceptanceTracker by \code{\link{adjustToAcceptanceRate}} before invoking this method
            popAcceptanceRates <- computePopulationAcceptanceRates(object@acceptanceTracker)
            nStepBackPops <- .getNStepBack(object, popAcceptanceRates)
            #
            nChainInPopulation <- getNChainInPopulation(sampleLogs)
            nChain <- nChainInPopulation*length(iPopulationsStep)
            zx <- array(NA_real_, dim=c(
                            length(object@iParametersThatRequireJumpProposal)
                            , nGeneration
                            , nChain
                            , 4L)
                        , dimnames=list(getParameterNames(sampleLogs)[object@iParametersThatRequireJumpProposal], NULL,NULL,NULL))
            acceptanceRatesChains <- numeric(nChain)    # acceptance rate of the chains group                
            fSampleStates <- if( isTRUE(object@is1DParameterVector) ) .sampleStates1D else .sampleStates
            for( iiPop in  seq_along(iPopulationsStep)){
                iPop <- iPopulationsStep[iiPop]
                # take care that iPop refers to Populations within all, iiPop only in iPopulationsStep
                # take care that iChainsForPopulationLocal refers to chains among iPopulationsStep only
                iChainsForPopulationLocal <- (iiPop-1L)*nChainInPopulation + (1L:nChainInPopulation)
                # when selecting states and calculating differences (and checking that they are different), 
                # then regard only the parameters that require updates
                Z <- getParametersAndInitialForPopulation(sampleLogs, iPop)[object@iParametersThatRequireJumpProposal,, ,drop=FALSE]
                # Z(nParm x nState x nChain) 
                mz <- iCurrentSamples[iPop] + getNInitialSamplePops(sampleLogs)[iPop]
                ##details<<
                ## The random samples are drawn within subgroups among all chains to support localization.
                #zxPop <- array(NA_real_, dim=c(dim(Z)[1], nGeneration, dim(Z)[3], 4L), dimnames=c(dimnames(Z)[1], NULL,NULL,NULL))
                chainsGroups <- .classifyChains(Z[,max(1,mz-nStepBackPops[iPop]):mz, ,drop=FALSE]
                    ,nGroups = getNChainGroupsInPopulation(sampleLogs)
                    ,nChainInPopulation = dim(Z)[3]
                )
                # note that chainsGroups and chainsGroup index chains within the current population
                for( iGroup in seq_along(chainsGroups) ){
                    chainsGroup <- chainsGroups[[iGroup]]
                    acceptanceRateGroup <- computeChainGroupAcceptanceRate(object@acceptanceTracker, iPop, chainsGroup)
                    acceptanceRatesChains[ iChainsForPopulationLocal[chainsGroup] ] <- acceptanceRateGroup 
                    nStepBack <- .getNStepBack(object, acceptanceRateGroup)
                    ZGroup <- Z[,,chainsGroup ,drop=FALSE]
                    zxGroup <- fSampleStates(Z=ZGroup,mZ=mz,nStepBack=nStepBack, nSample=nGeneration )
                    zx[,,iChainsForPopulationLocal[chainsGroup],] <- zxGroup
                }
            }
            # generate difference vectors and the rExtra
            genPropRes <- .generateXPropThin(object, zx)
            #length(iPopulationsStep)*mcSetup$nChainPop)
            # return entire parameter vector not only iParametersThatRequireJumpProposal, remaining parameters 0 in jump
            sDim <- getSampleDimensions(sampleLogs)
            # extend to array that includes parameters that do not require jump proposals
            jump=array( 0, dim=c(getNParameter(sDim),nGeneration,nChain), dimnames=list(getParameterNames(sDim),NULL,NULL))
            jump[ object@iParametersThatRequireJumpProposal,,] <- genPropRes$jump
            genPropRes$jump <- jump
            genPropRes
        })

 .classifyChains <- function(
         Z      ##<< numeric array (nParm, nSample, nChain): samples of the population
        , nGroups=1L
        , nChainInPopulation= dim(Z)[3]
){
    if( nGroups == 1L) return(list(group1=1:nChainInPopulation))
    # calculate eucledian distance between chains, with dimensions weighted by 1/variance inside chain
    currState <- adrop(Z[,ncol(Z), ,drop=FALSE],2L)
    #ZChain <- Z[,,4L]
    varDimsChains0 <- apply(Z,3, function(ZChain){
                apply(ZChain,1,var)
            })
    # constrain variance, to prevent high variabiliy chains to yield very loo distances
    varDimsChains0[ varDimsChains0==0 ] <- min(varDimsChains0[ varDimsChains0 > 0 ])
    varDimsChains <- array( pmin( quantile(varDimsChains0,0.6), varDimsChains0 ), dim=dim(varDimsChains0)) 
    fDist <- Vectorize(function(i,j){
        sum((currState[,j] - currState[,i])^2/varDimsChains[,i], na.rm=TRUE)
    })
    iVars <- 1:ncol(currState)
    dist <- outer(iVars, iVars, fDist)
    # distances from x-y may differ from y-x because of different variances
    # choose the minimum of these two
    #distS <- pmin(dist, t(dist))    
    #distS <- pmax(dist + t(dist))/2    
    distS <- (dist + t(dist))/2    
    # cluster the distance matrix with trying to achieve
    r1 <- 1
    fc <- suppressWarnings(fanny(as.dist(distS), k=nGroups, memb.exp=1+r1))
    while( (r1 < 4) && !fc$convergence["converged"] ){
        r1 <- r1*1.2
        fc <- suppressWarnings(fanny(as.dist(distS), k=nGroups, memb.exp=1+r1))
    }
    if( !fc$convergence["converged"] ) stop(".classifyChains: fanny classification did not converge")
    #fc$membership
    cl <- cl0 <- fc$clustering
    # for the goups with too few entries, select those from groups with more entries with highest probability
    cnt <- table(cl)
    minCnt <- min(cnt)
    iFew <- which(cnt == minCnt)
    nChainInGroup <- nChainInPopulation/nGroups
    iMore <- which(cnt > nChainInGroup)
    while( (minCnt < nChainInGroup) ){
         chainsOffer <- which(cl %in% iMore)
         offer <- fc$membership[chainsOffer, iFew, drop=FALSE ]
         posInd <- arrayInd( posScaler <- which.max(offer), dim(offer) )
         # if maximum probability of re-association is very low, break 
         prob <- offer[posInd]
         if( (min(cnt) > 1) && (prob < 0.1) ) break
         from <- chainsOffer[ posInd[1] ]
         to <- iFew[ posInd[2] ]
         cl[from] <- to
         cnt <- table(cl)
         minCnt <- min(cnt)
         iFew <- which(cnt == minCnt)
         iMore <- which(cnt > nChainInGroup)
    }
    #plot(currState[1,], currState[2,], type="n");  text( currState[1,], currState[2,], 1:ncol(currState), col=cl)
    ans <- lapply( 1:max(cl), function(i){which(cl==i)} )
    ans
 }
 
.sampleStates <- function(
        ### Sample records and initial state from past records.
        Z					##<< parameter sample (parms, steps,  chains)
        ,mZ					##<< number of steps in Z that have been sampled
        ,nStepBack		    ##<< integer scalar: number of recorded steps in recent past, used to generate the proposals
        ,nSample			##<< number of proposals to generate
#,rLogDen
){
    ##seealso<< 
    ## \code{\link{twDEMCBlockInt}}
    ## \code{\link{.generateXPropThin}}
    nChain = dim(Z)[3] 
    nParm = dim(Z)[1]
    ##details<<  
    ## Random states for chains for difference vector generation are only selecting within populations chains boundaries.
    ## This allows simulating independent population of chains.
    ## The acceptance rate may differ among populations. Hence, the set of previous generations
    ## , from which is randomly selected, may differ between poplations.
    # integer array (thinSteps*nChain*4) sample of chains, within populations
    #xC <- adrop(Z[,mZ,,drop=FALSE],2)		# assume current state xC as beginning of the interval: take last sample
    xC <- Z[,mZ,]		# assume current state xC as beginning of the interval: take last sample
    nStates <- nSample*nChain*3			# need to sample three states for snooker update
    seqChains <- 1:nChain 
    iStepsToSampleFrom <- max(1,(mZ-nStepBack+1)):mZ 
    rrStepsPop <-  sample(iStepsToSampleFrom, nStates, replace = TRUE)
    rrChainsPop <-  sample(seqChains, nStates, replace = TRUE)
    # in order to constrain two dimensions at the same time use the [] subset with an array see ?"["
    rr <- cbind(rrStepsPop,rrChainsPop)
    #zLogDen <- array(rLogDen[rr], dim=c(1,nSample*nChainPop,3), dimnames=list(parms="logDen", steps=NULL, zi=NULL) )
    rrParms <- cbind(parms=rep(1:nParm, nStates), steps=rep(rrStepsPop,each=nParm), chains=rep(rrChainsPop,each=nParm) )
    zParms <- array(Z[rrParms], dim=c(nParm,nSample*nChain,3), dimnames=list( rownames(Z), NULL, NULL) )
    # note that parameters are the first dimension		
    #rrParms <- cbind( rep(rrGenPop,nParm), rep(1:nParm, each=nStates), rep(rrChainsPop,nParm) )
    #zParms <- matrix( Z[rrParms], ncol=nParm)
    #z0 <- abind( zParms, zLogDen, along=1)
    #
    # append as forth initial state x vector to z
    # assume that chains is the last dimension, for each step assume the same initial state x 
    chainsZ <- rep(seqChains,nSample)
    Xs <- array(xC[,chainsZ], dim=c(nParm,nSample*nChain,1), dimnames=list(parms=rownames(Z), steps=NULL, zi=NULL) )
    #XsLogDen <- array( rLogDen[mZ,chainsZ],dim=c(1,nSample*nChainPop,1), dimnames=list(parms="logDen", steps=NULL, zi=NULL)  ) 
    #Xs <- abind( Xs, XsLogDen, along=1)	#logDen row never used for X
    #zx <- abind( zParms, Xs, along=3)
    zx <- array( c(zParms, Xs), dim=c(dim(zParms)[1:2], 4), dimnames=c(dimnames(zParms)[1], NULL,NULL) )
    #
    ##details<< 
    ## States z2 is  attempted to be distinct from z1.
    ## And z3 is attempted to be distinct from x (assuming steps in z are stacked chains collectives of x)
    iSame <- whichColsEqualSumHeuristics( zx[,,2], zx[,,1] )    # removed adrop for performance reasons, make sure other dimensions are not of size 1L
    nResample <- 8L
    i <- 0
    while( 0 < length(iSame) && i<nResample){
        zx[,iSame,2] <- zx[,sample.int(dim(zx)[2],length(iSame)),2]
        iSame <- whichColsEqualSumHeuristics( zx[,,2], zx[,,1] )    
        i<-i+1
    }
    if( !(i<nResample) && length(iSame)) warning(".sampleStates: unable to sample all distinct vectors for proposal generation.")
    #(tmp <- which( zx[,,2] == zx[,,1], arr.ind=TRUE))
    #when z3 is the same as x, i.e. zx[,,4]
    iSame <- whichColsEqualSumHeuristics( zx[,,3], zx[,,4] )
    #iSame <- unique(which( (zx[-nrow(zx),,3]==Xs), arr.ind=TRUE )[,2])
    i <- 0
    while( 0 < length(iSame) && i<nResample){
        zx[,iSame,3] <- zx[,sample.int(dim(zx)[2],length(iSame)),3]
        iSame <- whichColsEqualSumHeuristics( zx[,,3], zx[,,4] )
        i<-i+1
    }
    if( !(i<nResample) && length(iSame) ) warning(".sampleStates: unable to sample all distinct vectors for proposal generation.")
    #(tmp <- which( zx[,,3] == zx[,,4], arr.ind=TRUE))
    # reshape dimension of all steps into chains and steps
    # Second dimension so far are steps for all chains (anyway do not distinguish between chains of one population)
    # save reshaping, because it makes no difference to the assignment upstream in .sampleJumpsFromHistory
    zx
    #zxD <- array(zx, dim=c(nParm, nSample, nChain, 4L), dimnames=c(dimnames(zParms)[1], NULL,NULL,NULL))
    ##value<< random states (nParm, nSample, nChainPop, 4)
    ## First dimension is the state vector.
    ## Third (zx) dimension: three random vectors per step, forth, is the initial state x of the chain
    #zxD
}	

.sampleStates1D <- function(
        ### Sample records and initial state from past records.
        Z					##<< 1D parameter sample (steps,  chains)
        ,mZ					##<< number of steps in Z that have been sampled
        ,nStepBack		    ##<< integer scalar: number of recorded steps in recent past, used to generate the proposals
        ,nSample			##<< number of proposals to generate
#,rLogDen
){
    ##seealso<< 
    ## \code{\link{twDEMCBlockInt}}
    ## \code{\link{.generateXPropThin}}
    nChain = dim(Z)[3] 
    nParm = dim(Z)[1]
    ##details<<  
    ## Random states for chains for difference vector generation are only selecting within populations chains boundaries.
    ## This allows simulating independent population of chains.
    ## The acceptance rate may differ among populations. Hence, the set of previous generations to 
    ## randomly select from also may differ between poplations.
    # integer array (thinSteps*nChain*4) sample of chains, within populations
    #xC <- adrop(Z[,mZ,,drop=FALSE],2)		# assume current state xC as beginning of the interval: take last sample
    xC <- Z[,mZ,]
    nStates <- nSample*nChain*3			# need to sample three states for snooker update
    seqChains <- 1:nChain 
    iStepsToSampleFrom <- max(1,(mZ-nStepBack+1)):mZ 
    rrStepsPop <-  sample(iStepsToSampleFrom, nStates, replace = TRUE)
    rrChainsPop <-  sample(seqChains, nStates, replace = TRUE)
    # in order to constrain two dimensions at the same time use the [] subset with an array see ?"["
    rr <- cbind(rrStepsPop,rrChainsPop)
    #zLogDen <- array(rLogDen[rr], dim=c(1,nSample*nChainPop,3), dimnames=list(parms="logDen", steps=NULL, zi=NULL) )
    rrParms <- cbind(parms=rep(1:nParm, nStates), steps=rep(rrStepsPop,each=nParm), chains=rep(rrChainsPop,each=nParm) )
    #zParms <- array(Z[rrParms], dim=c(nParm,nSample*nChain,3), dimnames=list( parms=rownames(Z), steps=NULL, zi=NULL) )
    zParms <- array(Z[rrParms], dim=c(nSample*nChain,3) )
    # note that parameters are the first dimension		
    #rrParms <- cbind( rep(rrGenPop,nParm), rep(1:nParm, each=nStates), rep(rrChainsPop,nParm) )
    #zParms <- matrix( Z[rrParms], ncol=nParm)
    #z0 <- abind( zParms, zLogDen, along=1)
    #
    # append as forth initial state x vector to z
    # assume that chains is the last dimension, for each step assume the same initial state x 
    chainsZ <- rep(seqChains,nSample)
    Xs <- array(xC[chainsZ], dim=c(nSample*nChain,1), dimnames=list(steps=NULL, zi=NULL) )
    #XsLogDen <- array( rLogDen[mZ,chainsZ],dim=c(1,nSample*nChainPop,1), dimnames=list(parms="logDen", steps=NULL, zi=NULL)  ) 
    #Xs <- abind( Xs, XsLogDen, along=1)	#logDen row never used for X
    zx <- cbind( zParms, Xs)
    #
    ##details<< 
    ## States z2 is  attempted to be distinct from z1.
    ## And z3 is attempted to be distinct from x (assuming steps in z are stacked chains collectives of x)
    iSame <- which( zx[,2] == zx[,1] )
    nResample <- 8L
    i <- 0
    while( 0 < length(iSame) && i<nResample){
        zx[iSame,2] <- zx[sample.int(dim(zx)[2],length(iSame)),2]
        iSame <- which( zx[,2] == zx[,1] )
        i<-i+1
    }
    if( !(i<nResample) && length(iSame)) warning(".sampleStates: unable to sample all distinct vectors for proposal generation.")
    #(tmp <- which( zx[,,2] == zx[,,1], arr.ind=TRUE))
    #when z3 is the same as x, i.e. zx[,,4]
    iSame <- which( zx[,3] == zx[,4] )
    #iSame <- unique(which( (zx[-nrow(zx),,3]==Xs), arr.ind=TRUE )[,2])
    i <- 0
    while( 0 < length(iSame) && i<nResample){
        zx[iSame,3] <- zx[sample.int(dim(zx)[2],length(iSame)),3]
        iSame <- which( zx[,3] == zx[,4] )
        i<-i+1
    }
    if( !(i<nResample) && length(iSame) ) warning(".sampleStates: unable to sample all distinct vectors for proposal generation.")
    #(tmp <- which( zx[,,3] == zx[,,4], arr.ind=TRUE))
    ##value<< random states (nParm,(steps*nChainPop), 4)
    ## First dimension is the state vector.
    ## Second dimension are steps for all chains (anyway do not distinguish between chains of one population)
    ## Third (zx) dimension: three random vectors per step, forth, is the initial state x of the chain
    array(zx, dim=c(1L,dim(zx)), dimnames=list(rownames(Z),NULL,NULL))
}	

if(!exists(".generateXPropThin")) setGeneric(".generateXPropThin", function(object, ...) standardGeneric(".generateXPropThin"))
setMethod(".generateXPropThin", signature=c("DEJumpProposer"),
        function( object
                ### Generate Proposals from given random states.
                ,zx		##<< output of \code{\link{.sampleStates}} (maybe stacked states)
        ){
            ##seealso<< 
            ## \code{\link{twDEMCBlockInt}}
            ## \code{\link{.sampleStates}}
            ## \code{\link{.xStepSnooker}}
            ## \code{\link{.xStepParallel}}
            #
            nParm <- nrow(zx)
            nSample <- ncol(zx)
            nChain <- dim(zx)[3]
            # view dimensions step and chain as one steps to generate
            nStep <- nSample*nChain
            zxFlat <- array(zx, dim=c(nParm, nStep, 4L), dimnames=list(rownames(zx),NULL,NULL))
            nChainInPopulation <- nChain/length(object@gammaAdjustmentPops)
            gammaAdjustment <- rep(object@gammaAdjustmentPops, each=nSample*nChainInPopulation) 
            z <- zxFlat[,,1:3,drop=FALSE]	#three random state vectors per step
            xC <- adrop(zxFlat[,,4,drop=FALSE],3)	#initial state vector for step
            #
            # calculate the difference vectors between states
            diffVecAndRExtra <- matrix( NA_real_, nrow=nParm+1, ncol=nStep, dimnames=c(list(c(rownames(z),"rExtra")),list(NULL)) )
            boSnooker <- runif(nStep) < object@deSettings["pSnooker"]      # boolean vector: with probability pSnooker, do a snooker update instead (see terBraak08) for step at index
            if( 0 < sum(boSnooker) ){
                diffVecAndRExtra[,boSnooker] <- .xStepSnooker(object, z[,boSnooker,,drop=FALSE],xC[,boSnooker,drop=FALSE])
            }
            if( 0 < sum(!boSnooker) )
                diffVecAndRExtra[,!boSnooker] <- .xStepParallel(object, z[,!boSnooker,,drop=FALSE], gammaAdjustment = gammaAdjustment[!boSnooker] )	
            #
            # within the matrix the second, i.e. the last, dimension is Nsteps*nChain 
            # we can make this explicit by reshaping the array by just changing the dim attribute (only works on the last dimension)
            xStepAndExtra <- array(diffVecAndRExtra, dim=c(nParm+1,nSample,nChain), dimnames=list(parms=rownames(diffVecAndRExtra),steps=NULL,chains=NULL) )
            # numeric array (Npar+1,Nchains,Nsteps): difference vectors in parameter space for steps and chains
            # last row is the extra LogDen associated with snooker update
            xStep <- xStepAndExtra[-nrow(xStepAndExtra),,,drop=FALSE]			#dito
            rExtra <- adrop(xStepAndExtra[nrow(xStepAndExtra),,,drop=FALSE],1)			#second dim (columns step within Thinning interval)
            list(jump=xStep, logSnookerDensityMultiplier=rExtra)
            ##value<< List with components \describe{
            ## \item{jump}{numeric array (nPar x nStep x nChain): difference vectors in parameter space}
            ## \item{logSnookerDensityMultiplier}{numeric matrix (Nsteps x NChain): some extra LogDensity from snooker update}}
            ## nStep is from dimenstion of argument zx (usually deSettings["thin"])
        })

        
#if(!exists(".xStepSnooker")) setGeneric(".xStepSnooker", function(object, ...) standardGeneric(".xStepSnooker"))
#setMethod(".xStepSnooker", signature=c("DEJumpProposer"),
.xStepSnooker <- cmpfun(        
        function( object
                ### Generates Snooker updates based on given random numbers.
                ,z ##<< numeric array (Nparms,(nsteps), 3) of random states, dimnames parms,steps,zi
                ,xC ##<< current state (Nparms,(nsteps)) corresponding to chain of second dimension in z 
        ){
            # DE-Snooker update
            ##seealso<< 
            ## \code{\link{.generateXPropThin}}
            
            nParm <- nrow(z)
            nState <- ncol(z)
            #gamma_snooker =1.7 * some multiplicator to reach all states
            gamma_snooker = runif(nState, min=1.2,max=2.2)
            res <- matrix( as.numeric(NA), nrow=nParm+1, ncol=nState, dimnames=c(list(c(rownames(z),"rExtra")),list(NULL)) )
            #need loop because of inner products
            for( i in 1:nState){
                x_z = xC[,i] - z[,i,3]
                D2 = max(1.0e-300, x_z %*% x_z)
                projdiff = ((z[,i,1] -z[,i,2]) %*% x_z)/D2  # inner_product of difference with x_z / squared norm x_z
                res[-(nParm+1),i] <- xStepChain <-  (gamma_snooker[i] * projdiff) * x_z
                xPropChain = xC[,i] + xStepChain
                x_z = xPropChain - z[,i,3]
                D2prop = max((x_z%*%x_z), 1.0e-30)
                res[(nParm+1),i] <- rExtra <- object@deSettings["Npar12"] * (log(D2prop) - log(D2))   # Npar12  =(nParm - 1)/2  # extra term in logr for accept - reject ratio
            }
            res
        }) # DE-Snooker update

if(!exists(".xStepParallel")) setGeneric(".xStepParallel", function(object, ...) standardGeneric(".xStepParallel"))
setMethod(".xStepParallel", signature=c("DEJumpProposer"),
        function( object
                ### DE-parallel direction update based on given random numbers.
                ,z ##<< numeric array (Nparms,(nsteps), 3) of random states, dimnames parms,steps,zi
                ,gammaAdjustment=rep(1, ncol(z))  ##<< numeric vector (nSteps) of adjustment of multiplicator gamma to control acceptance rate
        ){
            ##seealso<< 
            ## \code{\link{.generateXPropThin}}
            nParm <- nrow(z)
            nState <- ncol(z)
            dz <- if( object@is1DParameterVector ){
                # need to reshape, because parameter dimension lost
                array( z[,,1] - z[,,2], dim=dim(z)[1:2], dimnames=dimnames(z)[1:2]) #jump vector as the difference between two random states (from z2 towards z1)
            } else {
                z[,,1] - z[,,2]
            }
            # in pGamma1 of the cases use a different multiplicate factor close to 1 to jump between modes
            gamma_par <- matrix( object@deSettings["F1"], nrow=nParm, ncol=nState) 
            boGammasF1 <- (runif(nState) < object@deSettings["pGamma1"])
            # in other cases scale the difference to get variability
            gamma_par[,!boGammasF1] <- object@deSettings["F2"] *
                    gammaAdjustment[!boGammasF1] *
                    runif( nParm*sum(!boGammasF1), min=1-object@deSettings["epsMult"], max=1+object@deSettings["epsMult"])    # multiplicative error to be applied to the difference 	
            xStepChain <- if (object@deSettings["epsAdd"] == 0) {  # avoid generating normal random variates for performance reasons (see terBraak08)
                        gamma_par * dz 
                    } else {
                        gamma_par * dz  + rnorm(length(dz),0,object@deSettings["epsAdd"])
                    }
            xStepChainAndRExtra <- rbind(xStepChain,rExtra=0)
            ##value<< Numeric matrix (nParm+1, nStep): difference vectors in parameter space.
            ## In addition to the parameter vector, the rExtra=0 is returned to denote that it is no Snooker update
            xStepChainAndRExtra
        })

if(!exists("setDeSetting")) setGeneric("setDeSetting", function(object, property,...) standardGeneric("setDeSetting"))
#' @export
setMethod("setDeSetting", signature=c("DEJumpProposer", "character"),
        function(object, property, value){
               object@deSettings[property] <- value
               object
        })


#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("DEJumpProposer")
if(!exists("isBurnin")) setGeneric("isBurnin", function(object,...) standardGeneric("isBurnin"))

setMethod("isBurnin", signature="DEJumpProposer", function(object,...
        ### Getter method for slot isBurnin
        ) {object@isBurnin})
if(!exists("isBurnin<-")) setGeneric("isBurnin<-", function(object,value) standardGeneric("isBurnin<-"))
setReplaceMethod("isBurnin", signature=c("DEJumpProposer", "logical"), function(object, value
        ### Setter method for slot isBurnin
        ) {object@isBurnin <- value; object})

if(!exists("getDeSettings")) setGeneric("getDeSettings", function(object,...) standardGeneric("getDeSettings"))
#' @export
setMethod("getDeSettings", signature="DEJumpProposer", function(object,...
        ### Getter method for slot deSettings
        ) {object@deSettings})
if(!exists("deSettings<-")) setGeneric("deSettings<-", function(object,value) standardGeneric("deSettings<-"))
#' @export
setReplaceMethod("deSettings", signature=c("DEJumpProposer", "DESettings"), function(object, value
        ### Setter method for slot deSettings
        ) {object@deSettings <- value; object})



