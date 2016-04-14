#' @include IntermediateSpecification.R
#' @include IntermediateUpdater.R

setClass("IntermediateUpdaters"
        ,representation(
                intermediateUpdaters="list"
        )
)

newIntermediateUpdaters <- function(intermediateSpecifications, parameterNames){
    initializeIntermediateUpdaters( new("IntermediateUpdaters"), intermediateSpecifications, parameterNames)
}

if(!exists("initializeIntermediateUpdaters")) setGeneric("initializeIntermediateUpdaters", function(object, intermediateSpecifications, parameterNames) standardGeneric("initializeIntermediateUpdaters"))
setMethod("initializeIntermediateUpdaters", signature=c("IntermediateUpdaters","list", "character"), 
        function(object,intermediateSpecifications, parameterNames) {
            if( length(intermediateSpecifications) ){
                specNames <- names(intermediateSpecifications)
                .errorOnNonUniqueNames( specNames, "intermediateSpecifications list" )
                object@intermediateUpdaters <- structure(lapply( specNames, function(specName){
                            spec <- intermediateSpecifications[[specName]]
                            iUpdater <- initializeIntermediateUpdater( new("IntermediateUpdater")
                                    ,specName, spec, parameterNames
                            )
                        }),names=specNames)
                object@intermediateUpdaters <- lapply(object@intermediateUpdaters, function(iUpdater){
                            initializeRequisiteUpdaters( iUpdater, object )
                        })
            }
            object
        })

if(!exists("getIntermediateUpdater")) setGeneric("getIntermediateUpdater", function(object,...) standardGeneric("getIntermediateUpdater"))
setMethod("getIntermediateUpdater", signature="IntermediateUpdaters", 
        function(object,intermediateId) {
            result <- object@intermediateUpdaters[[intermediateId]]
            assert_that(is(result,"IntermediateUpdater"))
            return( result)
        })

#if(!exists("getIntermediateUpdatersForIds")) setGeneric("getIntermediateUpdatersForIds", function(object,...) standardGeneric("getIntermediateUpdatersForIds"))
#setMethod("getIntermediateUpdatersForIds", signature="IntermediateUpdaters", 
        getIntermediateUpdatersForIds <- function(object,intermediateIds) {
            iPositions <- match(intermediateIds, names(object@intermediateUpdaters))
            iMissing <- which(is.na(iPositions))
            if( length(iMissing) )
                stopDemac("Folling intermediateIds are not within IntermediateUpdaters: "
                        ,intermediateIds[iMissing])
            result <- object@intermediateUpdaters[iPositions]
            return( result)
        }
   # )

if(!exists("getDependentIntermediateUpdaters")) setGeneric("getDependentIntermediateUpdaters", function(object,...) standardGeneric("getDependentIntermediateUpdaters"))
setMethod("getDependentIntermediateUpdaters", signature="IntermediateUpdaters", 
        function(object,
                ### Recursively return those intermetdiateUpdaters that depend on parameter of given position
                iParameters ##<< integer vector: parameter positions
        ) {
            directUpdaters <- lapply(object@intermediateUpdaters,function(intUpdater){
                        if( any(getBlockIndicesIParametersUsed(intUpdater) %in% iParameters)) return(intUpdater)
                        NULL
                    })  
            directUpdaters <- directUpdaters[ !vapply(directUpdaters, is.null, logical(1)) ]
            indirectUpdaters <- getIntermediateUpdatersDependingOnIntermediates( object, names(directUpdaters) )
            return( c(directUpdaters, indirectUpdaters) )
        })

if(!exists("getIntermediateUpdatersDependingOnIntermediates")) setGeneric("getIntermediateUpdatersDependingOnIntermediates", function(object,...) standardGeneric("getIntermediateUpdatersDependingOnIntermediates"))
setMethod("getIntermediateUpdatersDependingOnIntermediates", signature="IntermediateUpdaters", 
        function(object,intermediateIds) {
            newIds <- foundIds <- intermediateIds
            idsToTest <- setdiff( names(object@intermediateUpdaters), foundIds )
            while( length(newIds) && length(idsToTest)){
                newIds <- idsToTest[ vapply( object@intermediateUpdaters[idsToTest], function(intUpdater){
                            any(getIntermediatesUsed(intUpdater) %in% newIds)
                        }, logical(1)) ]
                foundIds <- c(foundIds, newIds)
                idsToTest <- setdiff( names(object@intermediateUpdaters), foundIds ) 
            }
            intIds <- setdiff(foundIds, intermediateIds)  # omit the requisite one
            object@intermediateUpdaters[ intIds ]            
        })


#library(twDev)    # automatic generation of GSetter
#--- generateAndPrintS4GSetters("IntermediateUpdaters")

if(!exists("getNames")) setGeneric("getNames", function(object) standardGeneric("getNames"))
setMethod("getNames", signature="IntermediateUpdaters", function(object) {names(object@intermediateUpdaters)})





