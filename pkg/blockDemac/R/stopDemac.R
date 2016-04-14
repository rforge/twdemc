#' @export
stopDemac <- function(
        ### Extension of \code{\link{stopCustom}} to generates error of given subclass of "demacError"
        ...                      ##<< further arguments concatenated to a message
        , subClass = character(0) ##<< string  vector: subClasses of error, (with most specific first)
        , call = sys.call(-1)     ##<< frame where the error occured 
) {
    stopCustom( c("demacError",subClass), call=call, ...)
}
attr(stopDemac,"ex") <- function(){
    myLog <- function(x) {
        if (!is.numeric(x))  stopDemac("myLog() needs numeric input, but input was",x)
        if (any(x < 0)) stopDemac(subClass="invalidValue", "myLog() needs positive inputs, but input was",x)
        log(x)
    }
    tryCatch(
            myLog(-3)
            #,invalidValue = function(condition) "invalid value"
            ,demacError = function(condition) paste("subclass of demacError:", condition$message)  
            ,error = function(condition) paste("general error:", condition$message)  
    )    
}

#' @export
stopDemacConvergenceProblems <- function(
        ### Extension of \code{\link{stopCustom}} to generates error of given subclass of c("demacError","demacConvergenceError") to signal convergence problems
        ...                      ##<< further arguments concatenated to a message
        , subClass = character(0) ##<< string  vector: subClasses of error, (with most specific first)
        , call = sys.call(-1)     ##<< frame where the error occured 
) {
    stopCustom( c(subClass,"demacConvergenceError","demacError"), call=call, ...)
}

#' @export
stopDemacInvalidChainState <- function(
        ### Extension of \code{\link{stopCustom}} to generates error of given subclass of c("demacError","demacInvalidChainState") to signal problems on computing intermediates 
        ...                      ##<< further arguments concatenated to a message
        , subClass = character(0) ##<< string  vector: subClasses of error, (with most specific first)
        , call = sys.call(-1)     ##<< frame where the error occured 
) {
    stopCustom( c(subClass,"demacInvalidChainState","demacError"), call=call, ...)
}




helloFromBlockDemac <- function(){
    msg <- "Hello World from blockDemac" 
    print(msg)
    msg
}

