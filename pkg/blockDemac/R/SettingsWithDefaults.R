#' @export
setClass("SettingsWithDefaults",
        contains="VIRTUAL",
        representation(  
                ctrl = "list"   
        )   
)

setMethod("show", "SettingsWithDefaults",
        function(object){
            str(object@ctrl)
        })

setMethod("initialize", signature=c("SettingsWithDefaults"),
        function( .Object, ...){
            .Object <- callNextMethod(.Object)
            .Object <- initializeDefaults(.Object)
            dots <- list(...)
            for(argName in names(dots)) {
                .Object[argName] <- dots[[argName]]
            }
            .Object
        })

setMethod("[", signature="SettingsWithDefaults",
        function(x,i,j,drop){
            if( !(i %in% names(x@ctrl)) ) stop("Tried getting undefined control parameter:",i)
            return( x@ctrl[[i]])
        })

setReplaceMethod("[", signature="SettingsWithDefaults",
        function(x,i,j,value){
            if( !(i %in% names(x@ctrl)) ) stop("Tried setting undefined control parameter:",i)
            x@ctrl[[i]] <- value
            validObject(x)
            return(x)
        })

if(!exists("getNames")) setGeneric("getNames", function(object) standardGeneric("getNames"))
setMethod("getNames", signature="SettingsWithDefaults", 
        function(object) {
            names(object@ctrl)
        }
)


if(!exists("initializeDefaults")) setGeneric("initializeDefaults", function(object,...) standardGeneric("initializeDefaults"))
