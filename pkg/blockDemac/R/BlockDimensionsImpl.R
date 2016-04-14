#' @include BlockDimensions.R
#' @include BlockSpecification.R

setClass("BlockDimensionsImpl",
        contains="BlockDimensions",
        representation(
                iBlockLogDensityComponents = "integer"	    # index name of of block by logDensityComponent (position)
                ,iBlockParametersUpdated = "integer"	    # index of block by parameter (position)
                ,blockNames = "character"	                # names of blocks
                ,isProposalRequired = "logical"	            # boolean (length(iBlockParametersUpdated): TRUE if proposal is required
        )
)

# mind to update blockDimensions<- when adding fields


newBlockDimensions <- function(blockSpecifications, parameterNames){
    initializeByBlockSpecifications(new("BlockDimensionsImpl")
            ,blockSpecifications = blockSpecifications, parameterNames=parameterNames)
}


if(!exists("initializeByBlockSpecifications")) setGeneric("initializeByBlockSpecifications", function(object, ...) standardGeneric("initializeByBlockSpecifications"))
setMethod("initializeByBlockSpecifications", signature=c("BlockDimensionsImpl"), 
        function(object,blockSpecifications, parameterNames) {
            specNames <- names(blockSpecifications)
            .errorOnNonUniqueNames( specNames, "blockSpecification list" )
            blockUpdaters <- structure(lapply(blockSpecifications, getBlockUpdater), 
                    names=specNames)
            object@blockNames <- specNames  
            object <- .initializeParameters(object, blockSpecifications, parameterNames ) # depends on densityComponentNames
            object <- .initializeLogDensityComponents(object, blockUpdaters)
            object <- .initializeProposalRequired(object, blockUpdaters)
            object
        }
)

setGeneric(".initializeParameters", function(object, ...) standardGeneric(".initializeParameters"))
setMethod(".initializeParameters", signature=c("BlockDimensionsImpl"), 
        function(object,blockSpecifications, parameterNames) {
            specNames <- names(blockSpecifications)
            object@iBlockParametersUpdated <- structure(rep( NA_integer_, length(parameterNames)), names=parameterNames)
            for( iBlock in 1:length(blockSpecifications)){
                bs <- blockSpecifications[[iBlock]]
                # if no parameters specified, assume all parameters updated
                iParametersToUpdate <- if( !length(getParametersToUpdate(bs)) ) seq_along(parameterNames) else {
                            .matchParameterNames(getParametersToUpdate(bs), parameterNames
                                    ,fErrorMsgList=function(parNames){list("unknown parameters ",parNames
                                                ," updated by blocks ",specNames[iBlock])} )                                
                        } 
                iDouble <- which(!is.na(object@iBlockParametersUpdated[iParametersToUpdate]))
                if( length(iDouble) ) stopDemac("block ",specNames[iBlock]
                            ," wants to update parameters ",parameterNames[iDouble]
                            ,", that are already updated by blocks ", specNames[object@iBlockParametersUpdated[iParametersToUpdate][iDouble]])
                object@iBlockParametersUpdated[iParametersToUpdate] <- iBlock                
            }
            iMissingParameters <- which(is.na(object@iBlockParametersUpdated))            
            if( length(iMissingParameters) ) stopDemac("Following parameters are not assigned to an updater: ", parameterNames[iMissingParameters] )
            object
        }
)


setGeneric(".initializeLogDensityComponents", function(object,...) standardGeneric(".initializeLogDensityComponents"))
setMethod(".initializeLogDensityComponents", signature=c("BlockDimensionsImpl"), 
        function(object, blockUpdaters) {
            logDensityComponentNamesBlocks <- lapply(blockUpdaters,function(blockUpdater){
                        getLogDensityComponentNames(blockUpdater)
                    })
            logDensityComponentNames <- as.vector(do.call( c,logDensityComponentNamesBlocks)) 
            nLogDensityComponents <- sapply(logDensityComponentNamesBlocks,length)
            object@iBlockLogDensityComponents <- rep( 1:length(blockUpdaters), times=nLogDensityComponents )
            names(object@iBlockLogDensityComponents) <- logDensityComponentNames 
            object
        }
)

setGeneric(".initializeProposalRequired", function(object,...) standardGeneric(".initializeProposalRequired"))
setMethod(".initializeProposalRequired", signature=c("BlockDimensionsImpl"), 
        function(object, blockUpdaters) {
            isProposalRequiredByBlock <- vapply(blockUpdaters, isJumpProposalRequired, logical(1) )
            object@isProposalRequired <- structure( isProposalRequiredByBlock[ object@iBlockParametersUpdated ]
                ,names=names(object@iBlockParametersUpdated) )
            object
        }
)

if(!exists("blockDimensions<-")) setGeneric("blockDimensions<-", function(object,value) standardGeneric("blockDimensions<-"))
setReplaceMethod("blockDimensions", signature=c("BlockDimensionsImpl", "BlockDimensions"), 
        function(object,value) {
            object@iBlockLogDensityComponents <- integer(getNLogDensityComponent(value))
            names(object@iBlockLogDensityComponents) <- getLogDensityComponentNames(value) 
            object@iBlockParametersUpdated <- integer(getNParameter(value))
            names(object@iBlockParametersUpdated) <-getParameterNames(value)
            nBlock <- getNBlock(value)
            object@blockNames <- getBlockNames(value)
            for( iBlock in 1:nBlock){
                object@iBlockLogDensityComponents[ getILogDensityComponentsByBlock(value,iBlock)] <- iBlock
                object@iBlockParametersUpdated[ getIParametersUpdatedByBlock(value,iBlock)] <- iBlock
            }
            object@isProposalRequired <- isProposalRequiredForParameters(value)
            object
        })

setMethod("getNBlock", signature="BlockDimensionsImpl", function(object) {length(object@blockNames)})
setMethod("getBlockNames", signature="BlockDimensionsImpl", function(object) {object@blockNames})

setMethod("getNParameter", signature="BlockDimensionsImpl", function(object) {length(object@iBlockParametersUpdated)})
setMethod("getNParameterWithProposal", signature="BlockDimensionsImpl", function(object) {sum(object@isProposalRequired)})
setMethod("getParameterNames", signature="BlockDimensionsImpl", function(object) {names(object@iBlockParametersUpdated)})
setMethod("getIParametersUpdatedByBlock", signature=c("BlockDimensionsImpl","integer"), function(object, blockIndex) {
            # due to selecting components by name, names of iBlockParametersUpdated are blockNames
            # expect to be names of parameters
            iParameters <- which(object@iBlockParametersUpdated==blockIndex)
            structure( iParameters, names=getParameterNames(object)[iParameters] )
        })
setMethod("getIParametersUpdatedByBlock", signature=c("BlockDimensionsImpl","character"), function(object, blockIndex) {
            iBlock <- which( blockIndex == getBlockNames(object))
            getIParametersUpdatedByBlock(object, iBlock)
        })
setMethod("isProposalRequiredForParameters", signature=c("BlockDimensionsImpl","integer"), function(object, parameterIndices) {
            object@isProposalRequired[parameterIndices]
        })
setMethod("isProposalRequiredForParameters", signature=c("BlockDimensionsImpl","character"), function(object, parameterIndices) {
            object@isProposalRequired[parameterIndices]
        })
setMethod("isProposalRequiredForParameters", signature=c("BlockDimensionsImpl","missing"), function(object, parameterIndices) {
            object@isProposalRequired
        })

setMethod("getNLogDensityComponent", signature="BlockDimensionsImpl", function(object) {length(object@iBlockLogDensityComponents)})
setMethod("getLogDensityComponentNames", signature="BlockDimensionsImpl", function(object) {names(object@iBlockLogDensityComponents)})
setMethod("getILogDensityComponentsByBlock", signature=c("BlockDimensionsImpl","integer"), function(object, blockIndex) {
            which(object@iBlockLogDensityComponents==blockIndex)
        })
setMethod("getILogDensityComponentsByBlock", signature=c("BlockDimensionsImpl","character"), function(object, blockIndex) {
            iBlock <- which( blockIndex == getBlockNames(object))
            getILogDensityComponentsByBlock(object, iBlock)
        })





