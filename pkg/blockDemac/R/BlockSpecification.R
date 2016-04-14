#' @include BlockUpdater.R

#' @export
blockSpec <- function(
        ### create a block specification entry
        parametersToUpdate  ##<< character vector of parameters that are updated by this block
        , parametersUsed=parametersToUpdate    ##<< character vector of parameters that are used by this block
        , blockUpdater      ##<< object of class BlockUpdater that handles the updating.
        , intermediatesUsed=character(0) ##<< name of the intermediatesUsed, used by this block 
## See e.g. \code{\link{MetropolisBlockUpdater}}
){
    ##details<<
    ## Omitting parametersToUpdate or specifying character(0), indicates updating of all used parameters
    ## Omitting parametersUsed indicates using only parametersToUpdate
    if( missing(parametersToUpdate) || !length(parametersToUpdate) ) 
        parametersToUpdate <- if( !missing(parametersUsed)) parametersUsed else character(0)
    ##value<< object of class BlockSpecification
    new("BlockSpecification"
            , parametersToUpdate=parametersToUpdate
            , parametersUsed=parametersUsed
            , blockUpdater=blockUpdater
            , intermediatesUsed=intermediatesUsed
    )
}

#' @export
setClass("BlockSpecification",
    representation(
        parametersToUpdate = "character"
        ,parametersUsed = "character"
        ,blockUpdater="BlockUpdater"
        ,intermediatesUsed="character"
    )
    ,validity = function(object){
        if( 
            length(object@parametersUsed) &&
            (!all(object@parametersToUpdate %in% object@parametersUsed)) 
        ) warning("Some parameters to update are not among used parameters.")
        return(TRUE)
    }
)


#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("BlockSpecification")
if(!exists("getParametersToUpdate")) setGeneric("getParametersToUpdate", function(object) standardGeneric("getParametersToUpdate"))
#' @export
setMethod("getParametersToUpdate", signature="BlockSpecification", function(object) {object@parametersToUpdate})

if(!exists("getParametersUsed")) setGeneric("getParametersUsed", function(object) standardGeneric("getParametersUsed"))
#' @export
setMethod("getParametersUsed", signature="BlockSpecification", function(object) {object@parametersUsed})

if(!exists("getBlockUpdater")) setGeneric("getBlockUpdater", function(object) standardGeneric("getBlockUpdater"))
#' @export
setMethod("getBlockUpdater", signature="BlockSpecification", function(object) {object@blockUpdater})
    
if(!exists("getIntermediatesUsed")) setGeneric("getIntermediatesUsed", function(object) standardGeneric("getIntermediatesUsed"))
#' @export
setMethod("getIntermediatesUsed", signature="BlockSpecification", function(object) {object@intermediatesUsed})
