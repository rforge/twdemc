#' @export
intermediateSpec <- function(
        ### create a block specification entry
        updateFunction      ##<< function to update the intermediate result with first argument the parameter vector
        , argsUpdateFunction=list()        ##<< list with arguments to updateFunction in addition to the parameters
        ,parameters=character(0)       ##<< character vector of parameter names that are used by this intermediate
        , intermediates=character(0)   ##<< character vector of names of intermediates that this intermediate depends on
){
    ##value<< object of class BlockSpecification
    new("IntermediateSpecification"
            , parametersUsed = parameters
            , intermediatesUsed = intermediates
            , updateFunction = updateFunction
            , argsUpdateFunction = argsUpdateFunction
    )
}

#' @export
setClass("IntermediateSpecification",
    representation(
        parametersUsed = "character"
        , intermediatesUsed = "character"	# the different name (..Ids.. missing) is important to not double slot name in DependenceIndices
        ,updateFunction="function"
        , argsUpdateFunction = "list"
    )
)

if(!exists("copy")) setGeneric("copy", function(object, source) standardGeneric("copy"))
setMethod("copy", signature=c("IntermediateSpecification","IntermediateSpecification"), 
        function(object,source) {
            object@parametersUsed <- getParametersUsed(source)
            object@intermediatesUsed <- getIntermediatesUsed(source)
            object@updateFunction <- getUpdateFunction(source)
            object@argsUpdateFunction <- getArgsUpdateFunction(source)
            object
        })

#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("IntermediateSpecification")

if(!exists("getParametersUsed")) setGeneric("getParametersUsed", function(object) standardGeneric("getParametersUsed"))
#' @export
setMethod("getParametersUsed", signature="IntermediateSpecification", function(object) {object@parametersUsed})

if(!exists("getIntermediatesUsed")) setGeneric("getIntermediatesUsed", function(object) standardGeneric("getIntermediatesUsed"))
#' @export
setMethod("getIntermediatesUsed", signature="IntermediateSpecification", function(object) {object@intermediatesUsed})

if(!exists("getUpdateFunction")) setGeneric("getUpdateFunction", function(object) standardGeneric("getUpdateFunction"))
#' @export
setMethod("getUpdateFunction", signature="IntermediateSpecification", function(object) {object@updateFunction})

if(!exists("getArgsUpdateFunction")) setGeneric("getArgsUpdateFunction", function(object) standardGeneric("getArgsUpdateFunction"))
#' @export
setMethod("getArgsUpdateFunction", signature="IntermediateSpecification", function(object) {object@argsUpdateFunction})
