#' @include SettingsWithDefaults.R

#' @export
setClass("DESettings", contains="SettingsWithDefaults")

setMethod("initializeDefaults", signature="DESettings", 
        function(object
            ### Set default settings.
        ){
            ##value<< Settings object with the following entries set:
            object@ctrl <- within(object@ctrl,{        
                        F = 2.38 		##<< related to multiplicative error (F2=F/sqrt(2*Npar), see eps.mult, Default: 2.38
                        F1 = 1          ##<< controls jumping between modes together with pGamma1, defaults to 1 (multiplicative factor) 
                        pSnooker= 0.1	##<< probability of a snooker update (others parallel updates), Default 0.1
                        pGamma1 = 0.1	##<< probability of jumping to state of another chain (different modes), Default 0.1
                        # epsMult =0.2	##<< >0 gives d-dimensional variation around gamma. 
                        epsMult =0.05	##<< >0 provides variation around gamma. 
                        ## It adds scaled uncorrelated noise to the proposal. Default: 0.2 , 
                        ## Its advantage over eps.add is 
                        ## that its effect scales with the differences of vectors in the population whereas eps.add does not. 
                        ## If the variance of a dimensions is close to 0, eps.mult gives smaller changes. , 
                        ## A uniformly distributed error, i.e. F2*runif(1+-epsMult*prop) multiplied to difference vector 
                        ## from parallel update 
                        epsAdd = 0	   	##<< >0 is needed to ensure that all positions in the space can be reached. 
                        ## For targets without gaps, it can set small or even to 0. Default 0 , 
                        ## sd of normally distributed error added to proposal by parallel or snooker update. 
                        nGenAcceptAverage = 24 ##<< number of generations back over which the acceptance rate is computed Default 16
                        probUpDir=0.5 	##<< probability of direction between two states of increasing Density, Default 0.5 
                        ## , Increasing this during burin may accelerate convergence
                        initialAcceptanceRate=0.25	##<< numeric scalar initially assumed acceptance rate (before being tracked by real acceptance)  
                        ## Used to calculate the number of generations backwards to sample from.
                        ## If only one value is given, then it is replicated to the entire matrix
                        DRgamma=0		##<< factor for reducing step length [0..1) in delayed rejection step, 0 means no DR step, Default 0
                        minPCompAcceptTempDecr=0.15  ##<< if acceptance rate drops below minPCompAcceptTempDecr+0.02 times this level, employ delayed rejection (DR), Default 0.15
                        pIndStep = 1.5 ##<< assume state to be independent, after on average about those number of accepted steps, Default 1.5
                        nPastGen = 10  ##<< factor for determining the number of recent past states to sample during burnin. Default 10 
                        ## It is multiplied by the number of parameters. Past generations are calculated 
                        ## by deviding by the number of chains per population
                        #,useLoadBalancing=FALSE	##<< if set to TRUE, submits each step separaetely to a free node, else each node gets an entire chain process 
                        #,freeMasterNode=FALSE	##<< if set to TRUE, no job is submitted to first node, so that this node can dispatch jobs without waiting,, see \code{\link[twSnowfall]{sfFArgsApplyDep}}
                        returnIntermediates=TRUE  ##<< set to FALSE if intermediate result is large and transferring it between slaves slows down calculation
                        #,useMultiT = FALSE	##<< whether to downscale Temperature of result components during burnin
                        #moved to block,TFix = vector("numeric",0)		##<< named numeric vector: result components for which temperature shoudl not change
                        useConditionalProposal = FALSE ##<< if set to TRUE, poposal steps are generated conditional on the state of the other blocks (experimental)
                        #,controlOverfittingMinNObs = Inf   ##<< scalar positive integer: minimum number of observations where bias is controlled. Need at least about 20 observations to work. Default Inf says no bias correction. Replaced by dInfo$maxLogDen
                        #thin = 4   # used in different places outside DEJumpGenerator       ##<< thinning interval: number of generations before a sample is recorded: to foster independency of samples
                        isAdjustingStepLengthToAcceptanceRate = TRUE    ##<< set to FALSE to avoid adjustment of stepLength to acceptance rate                        
                        ##end<<
                    })
            return(object)
        })


if(!exists("initializeDimensionalSettings")) setGeneric("initializeDimensionalSettings", function(object,...) standardGeneric("initializeDimensionalSettings"))
setMethod("initializeDimensionalSettings", signature="DESettings", 
        function(object
            ### initialize the constants that depende on number of parameters, chains, and thinning interval
            , nParameterWithProposal    ##<< scalar integer: number of parameters that require a step proposal
            , nChainInChainGroup        ##<< scalar integer: number of chains within on sampling group of chains (see details)
            , thin                      ##<< thinning interval, i.e. distance of steps, when a sample is recorded 
        ){
            ##details<< 
            ## Proposals are generated by differences among past states of a group of chains.
            ## The chains in one population can be rearranged into changing groups, 
            ## so that proposals are localized.
            if( missing(nParameterWithProposal) ) stop("Missing argument nParameter")
            if( nParameterWithProposal == 0) warning("DESettings: specified no parameters that require jump proposal.")
            if( missing(nChainInChainGroup) ) stop("Missing argument nChainInChainGroup")
            if( missing(thin) ) stop("Missing argument thin")
            nReprRows <- computeNRepresentativeRows(nParameterWithProposal, nChainInChainGroup)
            if( !is.finite(nReprRows) && (nParameterWithProposal != 0)) stop("nReprRows not finite")            
            object@ctrl <- within(object@ctrl,{
                        nRepresentativeRows <- nReprRows 
                        Npar12  <-(nParameterWithProposal - 1)/2   
                        F2 = object["F"]/sqrt(2*nParameterWithProposal)
                        nIntervalPerRange = round(nReprRows/2)  ##<< number of thinning intervals between recalculation of jumps
                        widthAcceptWindow = round(  object["nGenAcceptAverage"] / thin )
                    })
            return(object)
        })


#' @export
computeNRepresentativeRows <- function(
        ### Compute the number of rows necessary to to represent the n-dimensional parameter space
        nParameter,		        ##<< the number of parameters to estimate
        nChainInPopulation 	    ##<< the number of chains per population 
){
    ##seealso<<   
    ## \code{\link{initZtwDEMCNormal}}
    
    ##details<< see terBraak 2006 and 2008
    res <- as.integer(max(6L,ceiling((8L*nParameter)/nChainInPopulation)))
    res
    ### length of each chain so that each population is initialized with at least 8*nPar cases 
}
