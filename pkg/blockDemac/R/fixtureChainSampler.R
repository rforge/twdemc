
.replicateVector <- function(x, nColumn, colFac){
    if( missing(colFac) )
        return( matrix(x, nrow=length(x), ncol=nColumn, dimnames=list(names(x),NULL)) )
    if( length(colFac)==1) colFac <- rep(colFac,nColumn)    
    sapply( 1:nColumn, function(iCol){ x*colFac[iCol]})
}
attr(.replicateVector,"ex") <- function(){
    x <- c(a=1,b=2,c=3,d=4,e=5)
    .replicateVector(x, 3)
    .replicateVector(x, 3, 1:3)
}

# default arguments are designed to test Metropolis updater
# depening on (unchainging due to step=0) isM, m is updated
# with default step for m=7 and for c=9
.fixtureChainSampler <- function(
        x0 = c(a=0, isM=1, m=0, isC=1, c=0)
        ,stepSingle = c(a=NA, isM=0, m=7, isC=0, c=9)   # isX stays at initial 1    # TODO: remove a updated by Functional updater
        ,logDensityComponents0 =c(m1=-100, m2=-100, parms=-200) 
        ,thin = 4L
        ,nInterval = 3L
){
    # as.list(environment()) records the current values of function arguments, match.call does not return the defaults
    # if x0Chains is a single vector, than generate x0 for other chains by multiplying x0 by iChainInPop
    fx <- within( as.list(environment()),{
#                nChain <- length(nSamplePopulations) * nChainInPopulation 
#                colFac <- rep(1:nChainInPopulation,nChainInPopulation)
#                if( isMissingChainStatesArg ){
#                    if(!is.matrix(x0Chains)) x0Chains <- .replicateVector(x0Chains,nChain,colFac)
#                    stopifnot( nChain == ncol(x0Chains) )
#                    # maybe only x0Chains is given as a matrix but nChain is missing
#                    if(!is.matrix(logDensityComponents0Chains)) logDensityComponents0Chains <- .replicateVector(logDensityComponents0Chains,nChain,colFac)
#                    stopifnot( nChain == ncol(logDensityComponents0Chains) )
#                    chainStates <- lapply( 1:nChain, function(i){
#                                new("ChainState", parameters=x0Chains[,i], logDensityComponents=logDensityComponents0Chains[,i])
#                            })  
#                }  
                
                fUpdateBlockAddTerm <- function(
                        x, iParametersToUpdate
                        , lowerParBounds, upperParBounds, intermediates
                        , term
                ){
                    list(
                            isUpdated = TRUE 
                            ,xC = x[ iParametersToUpdate ] + term
                            ,intermediate = NULL    # actually dont return NULL, better list(), here testing handled
                    )
                }
                
                fLogDenConst <- function(
                        x
                        ,logDensityComponents       ##<< numeric vector (nLDensityComponents): 
                        ## already known logDensitys. Only calculate the positions that are NA
                        ,...
                        ,acceptedLogDensityComponents
                        ,logEnv=NULL        # usually in caller: new.env(parent = emptyenv())
                ){
                    # if first parameter is 1 then accept, else not
                    if( is.environment(logEnv) ){
                        logEnv$log <- paste0(logEnv$log,",x=",catNamedVector(x,3)," logDenC=",catNamedVector(logDensityComponents))
                    }
                    ret <- if (x[1]==1) {
                                acceptedLogDensityComponents    # accepted  
                            } else {
                                Inf*acceptedLogDensityComponents    # not accepted
                            }
                    ret
                }

                blocks <- list(
                    add7 = blockSpec("a",, new("FunctionBasedBlockUpdater", 
                            fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=7))),
                    # important that isM is the first parameter used
                    met1 = blockSpec(c("isM","m"),c("isM","m","c"),new("MetropolisBlockUpdater",
                            fLogDensity=fLogDenConst,
                            argsFLogDensity=list(acceptedLogDensityComponents = c(m1=-10)),
                            logDensityComponentNames = c("m1"))), 
                    # important that isC is the first parameter used
                    met2 = blockSpec(c("isC","c"),c("isC","c","m"),new("MetropolisBlockUpdater",
                            fLogDensity=fLogDenConst,
                            argsFLogDensity=list(acceptedLogDensityComponents = c(m2=-10,parms=-2)),
                            logDensityComponentNames = c("m2","parms"))) 
                )
                blockDimensions <- newBlockDimensions(blocks, names(x0))
                
                chainState <- new("ChainState", parameters=x0, logDensityComponents=logDensityComponents0)
                
                nGen <- thin*nInterval
                stepInfoRange <- new("StepInfoRangeImpl", step=matrix(
                                rep(stepSingle, nGen), nrow=length(x0), ncol=nGen, dimnames=list(names(x0),NULL))
                        ,nLogDensityComponents=getNLogDensityComponent(blockDimensions))
                
#                blockUpdaters <- newBlockUpdaters(blocks, names(x0))
#                
#                chainSampler <- new("ChainSamplerImpl", blockUpdaters=blockUpdaters, thin=thin) 
#                
#                chainSampler <- setRangeSpecs(chainSampler,
#                        chainState=chainState, nInterval=nInterval, stepInfoRange=stepInfoRange)
#                
#                sampleDimensions <- initializeSampleDimensionsImpl( new("SampleDimensionsImpl")
#                                        , blockDimensions=blockDimensions
#                                        , thin=thin, nSamplePopulations=nSamplePopulations
#                                        , nChainInPopulation=nChainInPopulation)
#                
#                deSettings <- initializeDimensionalSettings(new("DESettings"), nParameter=nParameter, nChainInPopulation = nChainInPopulation, thin=thin)
#                jumpProposer <- new("DEJumpProposer", deSettings=deSettings)
            })
    fx
}
