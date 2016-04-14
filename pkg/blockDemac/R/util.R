.errorOnNonUniqueNames <- function(specNames, listName){
    if( !length(specNames) || !all(nzchar(specNames)) ) stopDemac(
                "Need to provide unique names with ",listName,", but some components are missing names.")
    countsName <- table(specNames)
    if( any(countsName != 1) ) stopDemac("Need to provide unique names with ",listName,", but following names have duplicates: ",paste(names(countsName)[countsName != 1],collapse=","))
    invisible(NULL)
}

.matchParameterNames <- function(
        ### get the indices of given parameterNamesToMatch and stop if not found
        parameterNamesToMatch, 
        parameterNamesAll,
        fErrorMsgList = function(parameterNames){list("following parameterNames not found: ",parameterNames )}
){
    iParametersUsed <- match(parameterNamesToMatch, parameterNamesAll) 
    iMissing <- which(is.na(iParametersUsed))
    if(length(iMissing)) do.call( stopDemac, fErrorMsgList(parameterNamesToMatch[iMissing]))
    iParametersUsed
}

#' @export
unlistOrigNames <- function (x) 
{
    subNames <- as.vector(sapply(x,names))
    ans <- structure( unlist(x, recursive=FALSE, use.names=FALSE), names=subNames )
    ans
}

#from twMisc:
.isParallel <- function(
        ### isClusterRunning and number of processors != 1 
        cl=NULL
) {
    ##seealso<< \code{\link{parLapplySave}}, \code{\link{isClusterRunning}}, \code{\link{twMisc}}
    .isClusterRunning(cl) && length(cl) != 1
}
.isClusterRunning <- function (
        ### test if given cluster is running
        cl  ##<< cluster to test
){
    ##seealso<< \code{\link{parLapplySave}}, \code{\link{isParallel}}, \code{\link{twMisc}}
    ##value<< TRUE if clusterEval worked, FALSE if not
    tryCatch(any(unlist(parallel::clusterEvalQ(cl, TRUE))), error = function(err) { FALSE })
}






