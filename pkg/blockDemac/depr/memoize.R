
.memoize1 <- function(
        ### single slot cache for results of eval(block) indexed by parms     
        cacheEnv    ##<< the environment used as a cache (\code{new.env( hash=FALSE, parent=emptyenv() )})
        , parms     ##<< the key (compared using identical)
        , block     ##<< the code block to be evaluated
){
    if( is.environment(cacheEnv) ){
        if( identical( cacheEnv$parms, parms) ){
            #print("using cached")
            return( cacheEnv$blockResult )
        } else {
            cacheEnv$blockResult <- blockRes <- eval.parent(block )
            cacheEnv$parms <- parms
            return( blockRes )
        }
    } else eval.parent(block )
}
attr(.memoize1,"ex") <- function(){
    parms <- c(a=1, b=2)
    c <- 4
    cacheEnv <- new.env( hash=FALSE, parent=emptyenv() )
    parms <- c(a=1, b=2)
    (tmp2 <- .memoize1( cacheEnv, parms, {tmp = c*parms["a"]; tmp*parms["b"]} ))
    parms <- c(a=1, b=2)
    (tmp3 <- .memoize1( cacheEnv, parms, {tmp = c*parms["a"]; tmp*parms["b"]} ))
    parms <- c(a=2, b=2)
    (tmp3b <- .memoize1( cacheEnv, parms, {tmp = c*parms["a"]; tmp*parms["b"]} ))
    #
    cacheEnv <- NULL
    parms <- c(a=2, b=2)
    (tmp3b <- .memoize1( cacheEnv, parms, {tmp = c*parms["a"]; tmp*parms["b"]} ))
}

.memoizeN <- function(
        ### several slot cache for results of eval(block) indexed by parms     
        cacheEnv    ##<< the environment used as a cache (\code{new.env( hash=FALSE, parent=emptyenv() )})
        , parms     ##<< the key (compared using identical)
        , block     ##<< the code block to be evaluated
        , n = 16     ##<< the number of slots
){
    if( is.environment(cacheEnv) ){
        i = 1
        while( (i <= n) && !identical( cacheEnv$parms[[i]], parms) )  i = i+1
        if( i <= n ){
            #print(paste('using result of slot',i))
            return( cacheEnv$blockResult[[i]] )
        } else {
            blockRes <- eval.parent(block )
            if( !is.list(cacheEnv$parms) ){
                cacheEnv$parms <- list()
                cacheEnv$parms[[1]] <- parms
                cacheEnv$parms[n+1] <- NULL     # make n entries so that parms[[n]] is correct index
                cacheEnv$blockResult <- list()
                cacheEnv$blockResult[[1]] <- blockRes
            }else{
                i <- sample.int(n,1)
                cacheEnv$parms[[i]] <- parms
                cacheEnv$blockResult[[i]] <- blockRes
            }
            return( blockRes )
        }
    } else eval.parent(block )
}
attr(.memoizeN,"ex") <- function(){
    c <- 4
    #
    cacheEnv <- new.env( hash=FALSE, parent=emptyenv() )
    #undebug(.memoizeN)
    parms <- c(a=1, b=2)
    (tmp2 <- .memoizeN( cacheEnv, parms, {tmp <- c*parms["a"]; tmp*parms["b"]} ))
    parms <- c(a=2, b=2)
    (tmp2b <- .memoizeN( cacheEnv, parms2, {tmp <- c*parms["a"]; tmp*parms["b"]} ))
    #
    parms <- c(a=1, b=2)
    (tmp3 <- .memoizeN( cacheEnv, parms, {temp <- c*parms["a"]; tmp*parms["b"]} ))
    parms <- c(a=2, b=2)
    (tmp3b <- .memoizeN( cacheEnv, parms2, {temp <- c*parms["a"]; tmp*parms["b"]} ))
}

