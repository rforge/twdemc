#' @include ChainState.R
#' @include StepInfo.R
#' @include SubSpace.R
#' @include BlockDimensionsImpl.R
#' @include IntermediateUpdaters.R

# need to export to work on cluster
#' @export
setClass("BlockUpdaters"
        ,representation(
                blockDimensions="BlockDimensions"
                ,blockUpdaters="list"
                ,dependentBlockIndices="list"	    ##<< for each block, integer vector of block indices that depend on the parameters updated 
                ,dependentIntermediateIds="list"	##<< for each block, character vector of names of intermediate that depend on the parameters updated by the block
                ,requisiteIntermediateIds="list"	##<< for each block, character vector of names of intermediate that are required by the block
                ,blockIndexByParameter="integer"	##<< index of block for each index of parameter, ie. position in vector
                ,intermediateUpdaters="IntermediateUpdaters"
        )
        ,prototype(blockDimensions=new("BlockDimensionsImpl"))
)

newBlockUpdaters <- function(blockSpecifications, intermediateSpecifications, parameterNames){
    initializeBlockRelations( new("BlockUpdaters"), blockSpecifications, intermediateSpecifications, parameterNames)
}

if(!exists("initializeBlockRelations")) setGeneric("initializeBlockRelations", function(object, ...) standardGeneric("initializeBlockRelations"))
setMethod("initializeBlockRelations", signature=c("BlockUpdaters"), 
        function(object,blockSpecifications, intermediateSpecifications, parameterNames) {
            if( !length(blockSpecifications) ) stop("Need to specify at least one block")
            subSpace <- initializeIndexedBoundsForParameterNames(new("SubSpace"), parameterNames)
            object@blockUpdaters <- structure(lapply(blockSpecifications, function(blockSpec){
                                blockUpdater <- getBlockUpdater(blockSpec)
                                subSpace(blockUpdater) <- subSpace
                                blockUpdater
                            }), 
                    names=names(blockSpecifications))
            object@blockDimensions <- newBlockDimensions(blockSpecifications, parameterNames)
            object@intermediateUpdaters <- newIntermediateUpdaters(intermediateSpecifications, parameterNames)
            # order is important
            object <- .initializeBlockIndices(object, blockSpecifications, parameterNames ) # depends on densityComponentNames
            object <- .initializeBlockDependencies( object ) # depends on indices and component names
            object
        }
)


if(!exists(".initializeBlockIndices")) setGeneric(".initializeBlockIndices", function(object, ...) standardGeneric(".initializeBlockIndices"))
setMethod(".initializeBlockIndices", signature=c("BlockUpdaters"), 
        function(object,blockSpecifications, parameterNames) {
            bi <- new("BlockIndices")
            blockDim <- object@blockDimensions  # already initialized by blockSpecifications
            iParametersUsedBlocks <- .getIParametersUsedByBlock(object, blockSpecifications,parameterNames)
            object@blockIndexByParameter <- integer(getNParameter(blockDim))
            for( iBlock in 1:getLength(object)){
                bs <- blockSpecifications[[iBlock]]
                bi@iParametersToUpdate <- getIParametersUpdatedByBlock(blockDim, iBlock)
                object@blockIndexByParameter[bi@iParametersToUpdate] <- iBlock
                bi@iParametersUsed <-  iParametersUsedBlocks[[iBlock]]
                bi@iParametersToUpdateInBlock <-  match(bi@iParametersToUpdate, bi@iParametersUsed)
                bi@iLogDensityComponents <- getILogDensityComponentsByBlock(blockDim,iBlock)
                bi@intermediateIdsUsed <- getIntermediatesUsed(bs)
                blockIndices(object@blockUpdaters[[iBlock]]) <- bi
            }
            object
        }
)

if(!exists(".getIParametersUsedByBlock")) setGeneric(".getIParametersUsedByBlock", function(object, ...) standardGeneric(".getIParametersUsedByBlock"))
setMethod(".getIParametersUsedByBlock", signature=c("BlockUpdaters"),
        function(object,blockSpecifications, parameterNames) {
            iParametersUsedBlocks <- list()
            for( iBlock in 1:length(blockSpecifications)){
                bs <- blockSpecifications[[iBlock]]
                parameterNamesToMatch <- getParametersUsed(bs)
                iParametersUsed <- if( !length(parameterNamesToMatch) ){
                            ##details<< If no used parameters were specified for the block, 
                            ## then assume only those are used that are updated by the block.
                            getIParametersUpdatedByBlock(object@blockDimensions, iBlock)   
                        } else {
                            .matchParameterNames(parameterNamesToMatch, parameterNames)
                        } 
                iParametersUsedBlocks[names(blockSpecifications)[iBlock]] <- list(structure(iParametersUsed,names=parameterNamesToMatch))                
            }
            iParametersUsedBlocks
        }
)

if(!exists(".initializeBlockDependencies")) setGeneric(".initializeBlockDependencies", function(object,...) standardGeneric(".initializeBlockDependencies"))
setMethod(".initializeBlockDependencies", signature=c("BlockUpdaters"), 
        function(object) {
            object@dependentIntermediateIds <- 
                    lapply( 1:getLength(object), function(iBlock) {
                                blockUpdater <- getBlockUpdaterForIndex(object, iBlock)
                                names(.getDependentIntermediateUpdaters(object, blockUpdater))
                            })
            names(object@dependentIntermediateIds) <- getNames(object)
            object@requisiteIntermediateIds <- 
                    lapply( 1:getLength(object), function(iBlock) {
                                blockUpdater <- getBlockUpdaterForIndex(object, iBlock)
                                names(.getRequisiteIntermediateUpdaters(object, blockUpdater))
                            })
            names(object@requisiteIntermediateIds) <- getNames(object)
            #dependentNonRequisiteIntermediateIds <- lapply( seq_along(object@dependentIntermediateIds), function(iBlock){
            #                setdiff(object@dependentIntermediateIds[[iBlock]], object@requisiteIntermediateIds[[iBlock]])
            #            })
            object@dependentBlockIndices <- 
                    lapply( 1:getLength(object), function(iBlock) {
                                blockUpdater <- getBlockUpdaterForIndex(object, iBlock)
                                iParametersUpdated <- getBlockIndicesIParametersToUpdate(blockUpdater)
                                #twutz 1510: need to invalidate also blocks that have these intermediates as requisites
                                idsIntermediatesUpdated <- object@dependentIntermediateIds[[iBlock]] #dependentNonRequisiteIntermediateIds[[iBlock]]
                                .isDependentBlockIndex <- function(iDependentBlock){
                                    #twutz 1510: need to invalidate current block (set proposal is set before density is calculated)
                                    #if( iDependentBlock == iBlock) return(FALSE)
                                    dependentBlockUpdater <- getBlockUpdaterForIndex(object, iDependentBlock)       
                                    isDirectlyDependent <- any( getBlockIndicesIParametersUsed(dependentBlockUpdater) %in% iParametersUpdated  )
                                    isIndirectlyDependent <- any( getBlockIndicesIntermediateIdsUsed(dependentBlockUpdater) %in% idsIntermediatesUpdated )
                                    isDirectlyDependent || isIndirectlyDependent                             
                                }
                                which(vapply( 1:getLength(object), .isDependentBlockIndex, logical(1)))
                            })
            names(object@dependentBlockIndices) <- getNames(object)
            object
        }
)


if(!exists(".getDependentIntermediateUpdaters")) setGeneric(".getDependentIntermediateUpdaters", function(object,...) standardGeneric(".getDependentIntermediateUpdaters"))
setMethod(".getDependentIntermediateUpdaters", signature=c("BlockUpdaters"),
        function(object
                ### Get those intermediateUpdaters that are depending on the parameters updated by given blockUpdater
                , blockUpdater
        ) {
            iParameters <- getBlockIndicesIParametersToUpdate(blockUpdater)            
            getDependentIntermediateUpdaters(object@intermediateUpdaters,iParameters)
        })

if(!exists(".getRequisiteIntermediateUpdaters")) setGeneric(".getRequisiteIntermediateUpdaters", function(object,...) standardGeneric(".getRequisiteIntermediateUpdaters"))
setMethod(".getRequisiteIntermediateUpdaters", signature=c("BlockUpdaters"),
        function(object
                ### Recursively get those intermediateUpdaters that are required by the given blockUpdater
                , blockUpdater
        ) {
            # only relies on ids not on other dependencies
            result <- list()
            newRequisiteIds <- getBlockIndicesIntermediateIdsUsed(blockUpdater) # direct intermediate ids
            while(length(newRequisiteIds)){
                intUpdaters <- getIntermediateUpdatersForIds(object@intermediateUpdaters,newRequisiteIds)
                result <- c( result, intUpdaters )
                # add requisites of the current requisite
                newRequisiteIds <- do.call(c, lapply(intUpdaters, getIntermediatesUsed ))
                # remove requisites that have already been processed 
                newRequisiteIds <- setdiff(newRequisiteIds, names(result) )
            }
            result
        })

#if(!exists("computeUpdatedChainState")) setGeneric("computeUpdatedChainState",  function(chainState, object,stepInfo){ standardGeneric("computeUpdatedChainState") })
#setMethod("computeUpdatedChainState", signature=c("ChainState","BlockUpdaters","StepInfo"), 
computeUpdatedChainState <-  function( 
        ### One generation of updating all blocks for one chain
        chainState    ##<< current state of the chain
        ,object
        , stepInfo 
) {
    updatedChainState <- chainState
    for( i in 1:getLength(object) ){
        blockUpdater <- object@blockUpdaters[[i]]
        updatedChainState <- tmp <- .computeInvalidIntermediatesForBlock( updatedChainState, object, blockUpdater )
        # updateBlockInChainState needs to stay S4, because overleaded by blockUpdater
        # TODO reorder arguments: blockUpdater first
        updatedChainState <- updateBlockInChainState( updatedChainState, blockUpdater, object, stepInfo )
    }
    updatedChainState
}
#)

.tmp.f <- function(){
    tmp2 <- updateBlockInChainState( tmp, blockUpdater, object, stepInfo )
}


#if(!exists(".computeInvalidIntermediatesForBlock")) setGeneric(".computeInvalidIntermediatesForBlock",  function(chainState, object, blockUpdater){ standardGeneric(".computeInvalidIntermediatesForBlock") })
#setMethod(".computeInvalidIntermediatesForBlock", signature=c("ChainState","BlockUpdaters"), 
.computeInvalidIntermediatesForBlock <-  function( 
        ### One generation of updating all blocks for one chain
        chainState      ##<< current state of the chain
        ,object         ##<< blockUpdaters object that knows dependencies
        ,blockUpdater   ##<< the updater of the block to calculate intermediates for
) {
    intUpdaters <- getIntermediateUpdatersForIds(object@intermediateUpdaters, getBlockIndicesIntermediateIdsUsed(blockUpdater))
    for( intUpdater in intUpdaters){
        intId <- getIntermediatesIntermediateId(intUpdater)
        if( !length(getChainStatesIntermediate(chainState, intId)) ){
            chainState <- letUpdaterUpdateIntermediateInChainState( intUpdater, chainState, object@intermediateUpdaters@intermediateUpdaters )
        }
    }
    chainState
}
#)

if(!exists("computeInvalidIntermediates")) setGeneric("computeInvalidIntermediates",  function(object,...){ standardGeneric("computeInvalidIntermediates") })
setMethod("computeInvalidIntermediates", signature=c("BlockUpdaters"), 
        function( 
                ### Recompute all invalid intermediates in given chainState
                object
                ,chainState    ##<< current state of the chain
        ) {
            updatedChainState <- chainState
            for( i in 1:getLength(object) ){
                blockUpdater <- object@blockUpdaters[[i]]
                updatedChainState <- .computeInvalidIntermediatesForBlock( updatedChainState, object, blockUpdater )
            }
            updatedChainState
        }
)

if(!exists("computeInvalidChainStateDensities")) setGeneric("computeInvalidChainStateDensities",  function(object,chainState){ standardGeneric("computeInvalidChainStateDensities") })
setMethod("computeInvalidChainStateDensities", signature=c("BlockUpdaters","ChainState"), 
        ### Let the blockUpdaters compute the densities for all NA in logDensityComponents
        function( object, 
                chainState  ##<< current state of the chain 
        ) {
            for( i in 1:getLength(object) ){
                blockUpdater <- object@blockUpdaters[[i]]
                blockLogDensityComponents <- getBlockLogDensityComponents(chainState, blockUpdater)
                if( any(is.na(blockLogDensityComponents))){
                    chainState <- .computeInvalidIntermediatesForBlock( chainState, object, blockUpdater )
                    chainState <- computeChainStatesLogDensityComponents(chainState,blockUpdater)
                    # TODO remove checks when debugged
                    tmp <- getBlockLogDensityComponents(chainState, blockUpdater)
                    if( any(is.na(tmp)) ) stop("Need to calculate valid (non-NA) logDensities. Maybe set to -Inf")
                }
            }
            chainState
        }
)

if(!exists("invalidateDependents")) setGeneric("invalidateDependents",  function( object, chainState,  ... ){ standardGeneric("invalidateDependents") })
setMethod("invalidateDependents", signature=c("BlockUpdaters","ChainState"), 
        function( object
                , chainState    ##<< the chainState where dependencies should be invalidated
                ,  blockIndices ##<< which parameterblocks to invalidate
        ) {
            # usually called from ChainState$"blockParametersToUpdate<-"
            # replaced method dispatch for performance reasons
            #blockIndex <- getBlockIndexForParameters(object, getBlockIndicesIParametersToUpdate(blockIndices) )
            blockIndex <- object@blockIndexByParameter[ getBlockIndicesIParametersToUpdate(blockIndices)[1] ] 
            chainState <- .invalidateIntermediatesThatDependOn(chainState, object, blockIndex)
            chainState <- tmp <- .invalidateBlocksThatDependOn(chainState,object, blockIndex)
            chainState
        }
)

if(!exists("getBlockIndexForParameters")) setGeneric("getBlockIndexForParameters",  function( object,... ){ standardGeneric("getBlockIndexForParameters") })
setMethod("getBlockIndexForParameters", signature=c("BlockUpdaters"), 
        function( object
                ### return the blockIndex that is reponsible to update given parameters
                , iParametersToUpdate   ##<< index of parameters in parameter vector
        ) {
            object@blockIndexByParameter[ iParametersToUpdate[1] ]
        }
)



if(!exists(".invalidateBlocksThatDependOn")) setGeneric(".invalidateBlocksThatDependOn",  function( chainState, object, ... ){ standardGeneric(".invalidateBlocksThatDependOn") })
setMethod(".invalidateBlocksThatDependOn", signature=c("ChainState","BlockUpdaters"), 
        function( chainState, object,  blockIndex ) {
            for( dependentBlockUpdater in getDependentBlockUpdaters(object, blockIndex) ){
                chainState <- invalidateBlockLogDensities(chainState, dependentBlockUpdater)
            }
            chainState
        }
)

if(!exists("getDependentBlockUpdaters")) setGeneric("getDependentBlockUpdaters", function(object, blockIndex) standardGeneric("getDependentBlockUpdaters"))
setMethod("getDependentBlockUpdaters", signature=c("BlockUpdaters","integer"), 
        function(object, blockIndex) {
            lapply( object@dependentBlockIndices[[blockIndex]], function(dependentBlockIndex){
                        getBlockUpdaterForIndex(object, dependentBlockIndex)
                    })
        }
)
setMethod("getDependentBlockUpdaters", signature=c("BlockUpdaters","character"), 
        function(object, blockIndex) {
            lapply( object@dependentBlockIndices[[blockIndex]], function(dependentBlockIndex){
                        getBlockUpdaterForIndex(object, dependentBlockIndex)
                    })
        }
)

if(!exists(".invalidateIntermediatesThatDependOn")) setGeneric(".invalidateIntermediatesThatDependOn",  function( chainState, object, ... ){ standardGeneric(".invalidateIntermediatesThatDependOn") })
setMethod(".invalidateIntermediatesThatDependOn", signature=c("ChainState","BlockUpdaters"), 
        function( chainState, object,  blockIndex ) {
            intermediateIds <- object@dependentIntermediateIds[[blockIndex]]
            chainState <- .invalidateChainStatesIntermediates(chainState, intermediateIds)
            chainState
        }
)

if(!exists("getDependentIntermediateUpdaters")) setGeneric("getDependentIntermediateUpdaters", function(object, ...) standardGeneric("getDependentIntermediateUpdaters"))
setMethod("getDependentIntermediateUpdaters", signature=c("BlockUpdaters"), 
        function(object, blockIndex) {
            #assert_that( is.character(blockIndex) )
            intIds <- object@dependentIntermediateIds[[blockIndex]]
            getIntermediateUpdatersForIds(object@intermediateUpdaters, intIds)
        }
)


if(!exists("blockUpdaterApply")) setGeneric("blockUpdaterApply", function(object, FUN,...) standardGeneric("blockUpdaterApply"))
setMethod("blockUpdaterApply", signature=c("BlockUpdaters","function"), 
        function(object, FUN, simplify = TRUE, ...) {
            res <- sapply( 1:getLength(object), function(iBlock) {
                        blockUpdater <- object@blockUpdaters[[iBlock]] #getBlockUpdaterForIndex(object,iBlock)
                        FUN( blockUpdater, ... )            
                    }, simplify=simplify)
            names(res) <- getNames(object)
            res
        }
)

if(!exists("subSpace<-")) setGeneric("subSpace<-", function(object,value) standardGeneric("subSpace<-"))
setReplaceMethod("subSpace", signature=c("BlockUpdaters", "SubSpace"), 
        function(object, value) {
            for( i in 1:getLength(object) ){
                subSpace(object@blockUpdaters[[i]]) <- value
            }
            object
        })


if(!exists("getLength")) setGeneric("getLength", function(object) standardGeneric("getLength"))
setMethod("getLength", signature="BlockUpdaters", 
        function(object) {
            length(object@blockUpdaters)
        }
)

if(!exists("getNames")) setGeneric("getNames", function(object) standardGeneric("getNames"))
setMethod("getNames", signature="BlockUpdaters", 
        function(object) {
            names(object@blockUpdaters)
        }
)

if(!exists("getBlockUpdaterForName")) setGeneric("getBlockUpdaterForName", function(object, blockName) standardGeneric("getBlockUpdaterForName"))
setMethod("getBlockUpdaterForName", signature=c("BlockUpdaters","character"), 
        function(object, blockName) {
            if( length(blockName) != 1) stop("Must provide a single blockName in getBlockUpdaterForName(blockUpdaters)")
            blockUpdater <-  object@blockUpdaters[[blockName]]
            if(is.null(blockUpdater)) stop("BlockName '",blockName,"' not found in blockUpdaters.")
            blockUpdater
        }
)

if(!exists("getBlockUpdaterForIndex")) setGeneric("getBlockUpdaterForIndex", function(object, index) standardGeneric("getBlockUpdaterForIndex"))
setMethod("getBlockUpdaterForIndex", signature=c("BlockUpdaters","integer"), 
        function(object, index) {
            if( length(index) != 1) stop("Must provide a single index in getBlockUpdaterForIndex(blockUpdaters)")
            if( (index < 1) || (index > getLength(object)) ) stop("Index out of range: ",index)        
            object@blockUpdaters[[index]]
        }
)

if(!exists("getIParametersThatRequireJumpProposal")) setGeneric("getIParametersThatRequireJumpProposal", function(object) standardGeneric("getIParametersThatRequireJumpProposal"))
setMethod("getIParametersThatRequireJumpProposal", signature=c("BlockUpdaters"), 
        function(object) {
            if( !length(object@blockUpdaters) ) stop("BlockUpdaters must be set before accessing other properties.")
            iParametersBlocks <- lapply( object@blockUpdaters, function(blockUpdater){
                        if( isJumpProposalRequired(blockUpdater) ) return( getBlockIndicesIParametersToUpdate(blockUpdater))
                        integer(0)
                    })
            names(iParametersBlocks) <- NULL
            sort(do.call(c,iParametersBlocks))
        }
)

if(!exists("getIntermediateIdsUsed")) setGeneric("getIntermediateIdsUsed", function(object,...) standardGeneric("getIntermediateIdsUsed"))
setMethod("getIntermediateIdsUsed", signature=c("BlockUpdaters"), function(object) {getNames(object@intermediateUpdaters)})


#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("BlockUpdaters")
if(!exists("getLogDensityComponentNames")) setGeneric("getLogDensityComponentNames", function(object) standardGeneric("getLogDensityComponentNames"))
setMethod("getLogDensityComponentNames", signature="BlockUpdaters", function(object) {stop("getLogDensityComponentNames: use BlockDimensions insted of BlockUpdaters")})

if(!exists("getParameterNames")) setGeneric("getParameterNames", function(object) standardGeneric("getParameterNames"))
setMethod("getParameterNames", signature="BlockUpdaters", function(object) {stop("getParameterNames: use BlockDimensions insted of BlockUpdaters")})

if(!exists("getBlockDimensions")) setGeneric("getBlockDimensions", function(object) standardGeneric("getBlockDimensions"))
setMethod("getBlockDimensions", signature="BlockUpdaters", function(object) {object@blockDimensions})

if(!exists("getIntermediateUpdaters")) setGeneric("getIntermediateUpdaters", function(object) standardGeneric("getIntermediateUpdaters"))
setMethod("getIntermediateUpdaters", signature="BlockUpdaters", function(object) {object@intermediateUpdaters})

