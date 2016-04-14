#' @include DependenceIndices.R

#------------------------- public data class used by BlockUpdater to represent indices
#' @export
setClass("BlockIndices",        
        ### Indices of subRanges used by a block in overall parameter and logDensity vectors
        contains = c("DependenceIndices"),             
        representation(  
                iParametersToUpdate="integer"	##<< indices within parameter vector
                     ## which parametes are updated by this block
                ,iParametersToUpdateInBlock="integer"	##<< indices of parameters to update within block parameter vector (that also includes used parameters) 
                ,iLogDensityComponents="integer"   ##<< indices within vector of combined logDensity results
        )
)

if(!exists("copyBlockIndices")) setGeneric("copyBlockIndices", function(object, source) standardGeneric("copyBlockIndices"))
#' @export
setMethod("copyBlockIndices", signature=c("BlockIndices","BlockIndices")
        , function(object,source) {
            object <- copyDependenceIndices(object, source)
            object@iParametersToUpdate <- getBlockIndicesIParametersToUpdate(source)
            object@iParametersToUpdateInBlock <- getBlockIndicesIParametersToUpdateInBlock(source)
            object@iLogDensityComponents <- getBlockIndicesILogDensityComponents(source)
            object
        })


##' @export
#setMethod("getBlockIndicesIParametersToUpdate", signature="BlockIndices", function(object) {object@iParametersToUpdate})
##' @export
#setMethod("getIParametersToUpdateInBlock", signature="BlockIndices", function(object) {
#            # bottleneck (called very often, hence omit check for length)
##            if( length(object@iParametersToUpdateInBlock) != length(object@iParametersToUpdate)) 
##                stopDemac("iParametersToUpdateInBlock not set yet.")
#            object@iParametersToUpdateInBlock
#        })
##' @export
#setMethod("getBlockIndicesILogDensityComponents", signature="BlockIndices", function(object) {object@iLogDensityComponents})


setMethod("show", "BlockIndices",
        function(object){
            cat("iParametersToUpdate=",object@iParametersToUpdate,sep="")
            cat(", iParametersUsed=",object@iParametersUsed,sep="")
            cat(", iLogDensityComponents=",object@iLogDensityComponents,sep="")
            cat(", intermediateIdsUsed='",paste(object@intermediateIdsUsed,collapse=","),"'",sep="")
            cat("\n")
        })

# define shortcuts to avoid costly method dispatch in R
getBlockIndicesIParametersToUpdate <- function(object) { object@iParametersToUpdate }
getBlockIndicesIParametersToUpdateInBlock <- function(object) { object@iParametersToUpdateInBlock }
getBlockIndicesILogDensityComponents <- function(object) { object@iLogDensityComponents }
