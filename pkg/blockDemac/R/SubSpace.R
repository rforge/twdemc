#' @export
setClass("SubSpace"
        , representation(  
                lowerParBounds="numeric"
                ,upperParBounds="numeric"
                ,splits = "numeric"
                ,indexedParBounds="matrix"	 
        )
        ,prototype = list(
                lowerParBounds=numeric(), upperParBounds=numeric(),
                indexedParBounds = matrix(numeric(0)),
                splits= numeric())
)


if(!exists("getSplitSubSpaces")) setGeneric("getSplitSubSpaces", function(object,splitPoint,...) standardGeneric("getSplitSubSpaces"))
#' @export
setMethod("getSplitSubSpaces", signature=c("SubSpace","numeric"), function(object,splitPoint,...) {
            object@splits <- c(object@splits, splitPoint)
            lower <- upper <- object
            lower@upperParBounds <- c(lower@upperParBounds, splitPoint)
            upper@lowerParBounds <- c(upper@lowerParBounds, splitPoint)
            list(
                    lower = initializeIndexedBoundsForParameterNames(lower, rownames(object@indexedParBounds)),
                    upper = initializeIndexedBoundsForParameterNames(upper, rownames(object@indexedParBounds))
            )
        })

if(!exists("isInSubSpace")) setGeneric("isInSubSpace", function(object, parameters) standardGeneric("isInSubSpace"))
#' @export
setMethod("isInSubSpace", signature=c("SubSpace","numeric"), 
        function(object, parameters) {
            namesConstrainedParameters <- intersect( names(parameters), names(object@upperParBounds) )
            constrainedParameters <- parameters[ namesConstrainedParameters ]
            upperConstraints <- object@upperParBounds[namesConstrainedParameters ] 
            if (!all(constrainedParameters <=  upperConstraints)) return( FALSE)
            #
            namesConstrainedParameters <- intersect( names(parameters), names(object@lowerParBounds) )
            constrainedParameters <- parameters[ namesConstrainedParameters ]
            lowerConstraints <- object@lowerParBounds[namesConstrainedParameters ] 
            if (!all(constrainedParameters >=  lowerConstraints)) return( FALSE)
            #
            return(TRUE)
        }
)



# possibilit to call without S4
.letSubSpaceReflectBoundaries <- function(object, parameters) {
    namesConstrainedParameters <- intersect( names(parameters), names(object@upperParBounds) )
    constrainedParameters <- parameters[ namesConstrainedParameters ]
    upperConstraints <- object@upperParBounds[namesConstrainedParameters ]
    overshoot <- constrainedParameters - upperConstraints
    parameters[namesConstrainedParameters] <- parameters[namesConstrainedParameters] - 
            ifelse( is.finite(overshoot) & (overshoot > 0),  2*overshoot, 0 )
    #
    namesConstrainedParameters <- intersect( names(parameters), names(object@lowerParBounds) )
    constrainedParameters <- parameters[ namesConstrainedParameters ]
    lowerConstraints <- object@lowerParBounds[namesConstrainedParameters ] 
    overshoot <- lowerConstraints - constrainedParameters
    parameters[namesConstrainedParameters] <- parameters[namesConstrainedParameters] + 
            ifelse( is.finite(overshoot) & (overshoot > 0),  2*overshoot, 0 )
    #
    parameters
}

.letSubSpaceReflectBoundariesMatrix <- function(
        ### reflect parameters given as a matrix on parameter bounds 
        object          ##<< subspace object
        , parameters    ##<< numeric matrix (nParm x nCases): parameters, with parameter names in column names
) {
    if( !is.matrix(parameters) ) stop(".letSubSpaceReflectBoundariesMatrix: must provide a matrix of parameters with column names")
    namesConstrainedParameters <- intersect( colnames(parameters), names(object@upperParBounds) )
    constrainedParameters <- parameters[, namesConstrainedParameters ]
    upperConstraints <- matrix(object@upperParBounds[namesConstrainedParameters ]
            , byrow=TRUE, ncol=length(namesConstrainedParameters), nrow=nrow(parameters))
    overshoot <- constrainedParameters - upperConstraints
    parameters[,namesConstrainedParameters] <- parameters[,namesConstrainedParameters] - 
            ifelse( is.finite(overshoot) & (overshoot > 0),  2*overshoot, 0 )
    #
    namesConstrainedParameters <- intersect( colnames(parameters), names(object@lowerParBounds) )
    constrainedParameters <- parameters[, namesConstrainedParameters ]
    lowerConstraints <- matrix(object@lowerParBounds[namesConstrainedParameters ] 
            , byrow=TRUE, ncol=length(namesConstrainedParameters), nrow=nrow(parameters))
    overshoot <- lowerConstraints - constrainedParameters
    parameters[,namesConstrainedParameters] <- parameters[,namesConstrainedParameters] + 
            ifelse( is.finite(overshoot) & (overshoot > 0),  2*overshoot, 0 )
    #
    parameters
}


if(!exists("reflectBoundaries")) setGeneric("reflectBoundaries", function(object, parameters) standardGeneric("reflectBoundaries"))
#' @export
setMethod("reflectBoundaries", signature=c("SubSpace","numeric"), .letSubSpaceReflectBoundaries )

if(!exists("initializeIndexedBoundsForParameterNames")) setGeneric("initializeIndexedBoundsForParameterNames", function(object,...) standardGeneric("initializeIndexedBoundsForParameterNames"))
#' @export
setMethod("initializeIndexedBoundsForParameterNames", signature=c("SubSpace"), function(object, parameterNames) {
            ##details<< 
            ## Construct a matrix with rows corresponding to given parameter names.
            ## first column gives the lower par bounds (-Inf for no bound),
            ## second column gives the upper par bounds (+Inf for no bound).
            ## This will speed up successive comparisons of parameter vectors.
            ##seealso<< \code{\link{getIndexedParBounds}}, \code{\link{isInIndexedParBounds}}
            object@indexedParBounds <- matrix(-Inf, nrow=length(parameterNames), ncol=2, dimnames=list(parameterNames,c("lower","upper")))
            object@indexedParBounds[,2L] <- +Inf
            namesConstrained <- intersect( parameterNames, names(object@lowerParBounds) )
            object@indexedParBounds[namesConstrained,1L] <- object@lowerParBounds[namesConstrained]
            namesConstrained <- intersect( parameterNames, names(object@upperParBounds) )
            object@indexedParBounds[namesConstrained,2L] <- object@upperParBounds[namesConstrained]
            object
        })

if(!exists("isInIndexedParBounds")) setGeneric("isInIndexedParBounds", function(object, parameters, ...) standardGeneric("isInIndexedParBounds"))
#' @export
setMethod("isInIndexedParBounds", signature=c("SubSpace","numeric"), 
        function(object, parameters, iParameters=TRUE) {
            all( parameters >= object@indexedParBounds[iParameters,1L] ) &
                    all( parameters <= object@indexedParBounds[iParameters,2L] )
        }
)

if(!exists("reflectAtIndexedBounds")) setGeneric("reflectAtIndexedBounds", function(object, parameters, ...) standardGeneric("reflectAtIndexedBounds"))
#' @export
setMethod("reflectAtIndexedBounds", signature=c("SubSpace","numeric"),
        function(object, parameters, iParameters=TRUE){
            lowerConstraints <- object@indexedParBounds[iParameters,1L] 
            isLowerConstrained <- is.finite(lowerConstraints)
            overshoot <- lowerConstraints - parameters
            isOvershoot <- overshoot > 0
            parameters[isLowerConstrained][isOvershoot] <- parameters[isLowerConstrained][isOvershoot] + 2*overshoot[isOvershoot]
            #
            upperConstraints <- object@indexedParBounds[iParameters,2L]
            isUpperConstrained <- is.finite(upperConstraints)
            overshoot <- parameters - upperConstraints
            isOvershoot <- overshoot > 0
            parameters[isUpperConstrained][isOvershoot] <- parameters[isUpperConstrained][isOvershoot] - 2*overshoot[isOvershoot]
            #
            # if correction undershot lower bound, set hard
            parameters[isLowerConstrained] <- pmax(lowerConstraints, parameters[isLowerConstrained])
            parameters
        })


#library(twDev)    # automatic generation of GSetter
#------------ generateAndPrintS4GSetters("SubSpace")

# only getter methods, hide implementation
if(!exists("getLowerParBounds")) setGeneric("getLowerParBounds", function(object) standardGeneric("getLowerParBounds"))
#' @export
setMethod("getLowerParBounds", signature="SubSpace", function(object) {object@lowerParBounds})

if(!exists("getUpperParBounds")) setGeneric("getUpperParBounds", function(object) standardGeneric("getUpperParBounds"))
#' @export
setMethod("getUpperParBounds", signature="SubSpace", function(object) {object@upperParBounds})

if(!exists("getIndexedParBounds")) setGeneric("getIndexedParBounds", function(object) standardGeneric("getIndexedParBounds"))
#' @export
setMethod("getIndexedParBounds", signature="SubSpace", function(object) {object@indexedParBounds})

#if(!exists("getSplits")) setGeneric("getSplits", function(object) standardGeneric("getSplits"))
#setMethod("getSplits", signature="SubSpace", function(object) {object@splits})

