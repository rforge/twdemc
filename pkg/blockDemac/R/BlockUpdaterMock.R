#' @include BlockUpdater.R
#' @include BlockUpdaters.R

# Mock that does modify the chainstate
# but gives logDensityComponents - used for initializing a BlockDimensions object 

setClass("BlockUpdaterMock"
    , contains ="BlockUpdater"	    # blockIndices, subSpace, chainState
    , representation(  
        logDensityComponentNames="character"
    )
)


setMethod("show", "BlockUpdaterMock",
        function(object){
            callNextMethod()
            print(object@logDensityComponentNames)
            #cat("chainState: "); show(object@chainState)
        })



setMethod("updateBlockInChainState", signature=c("ChainState","BlockUpdaterMock","BlockUpdaters","StepInfo")
    , function(chainState, blockUpdater, blockUpdaters, stepInfo) {
        chainState
    })

setMethod("getLogDensityComponentNames", signature="BlockUpdaterMock",  function(object) { object@logDensityComponentNames })

setMethod("isJumpProposalRequired", signature="BlockUpdaterMock",
    # assume when logDensity is reported, then also jumps are required
    function(object) {
        (length(object@logDensityComponentNames) != 0)
    }
)



blockSpecMock <- function(parametersToUpdate, parametersUsed, logDensityComponentNames){
    ##details<<
    ## Omitting parametersToUpdate or specifying character(0), indicates updating of all used parameters 
    ## Omitting parametersUsed or specifying character(0), indicates usage of all parameters 
    if( missing(parametersUsed) ) parametersUsed <- character(0)
    if( missing(parametersToUpdate) || !length(parametersToUpdate) ) parametersToUpdate <- parametersUsed
    if( missing(logDensityComponentNames) || !length(logDensityComponentNames) ) logDensityComponentNames <- character(0)
    blockUpdaterMock <- new("BlockUpdaterMock", logDensityComponentNames=logDensityComponentNames)
    new("BlockSpecification", parametersToUpdate=parametersToUpdate
            , parametersUsed=parametersUsed, blockUpdater=blockUpdaterMock )
}


