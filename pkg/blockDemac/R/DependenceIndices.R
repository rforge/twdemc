setClass("DependenceIndices",
        representation(
                iParametersUsed="integer"	    ##<< indices within parameter vector
                ## which parameters are used by this dependent object
                ,intermediateIdsUsed="character"     ##<< keys to intermediates list that an object depends on  
        )
)

if(!exists("copyDependenceIndices")) setGeneric("copyDependenceIndices", function(object, source) standardGeneric("copyDependenceIndices"))
#' @export
setMethod("copyDependenceIndices", signature=c("DependenceIndices","DependenceIndices")
        , function(object,source) {
            object@iParametersUsed <- getBlockIndicesIParametersUsed(source)
            object@intermediateIdsUsed <- getBlockIndicesIntermediateIdsUsed(source)
            object
        })


# define shortcuts to avoid costly method dispatch in R
getBlockIndicesIParametersUsed <- function(object) { object@iParametersUsed }
getBlockIndicesIntermediateIdsUsed <- function(object) { object@intermediateIdsUsed }





