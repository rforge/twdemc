DLLfuncTest <- function (
        ### Adjustment of \code{\link[deSolve]{DLLfunc}}, that does not coerce parms to double vector.
        func, times, y, parms, dllname, initfunc = dllname, 
        rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL, 
        initforc = NULL, fcontrol = NULL) 
{
    if (!is.numeric(y)) 
        stop("`y' must be numeric")
    n <- length(y)
    if (!is.null(times) && !is.numeric(times)) 
        stop("`times' must be NULL or numeric")
    if (!is.null(outnames)) 
        if (length(outnames) != nout) 
            stop("length outnames should be = nout")
    ModelInit <- NULL
    Outinit <- NULL
    flist <- list(fmat = 0, tmat = 0, imat = 0, ModelForc = NULL)
    Ynames <- attr(y, "names")
    if (is.null(dllname) || !is.character(dllname)) 
        stop("`dllname' must be a name referring to a dll")
    if (!is.null(initfunc)) 
        if (is.loaded(initfunc, PACKAGE = dllname, type = "") || 
                is.loaded(initfunc, PACKAGE = dllname, type = "Fortran")) {
            ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
        }
        else if (initfunc != dllname && !is.null(initfunc)) 
            stop(paste("cannot integrate: initfunc not loaded ", 
                            initfunc))
    if (is.null(initfunc)) 
        initfunc <- NA
    #if (!is.null(forcings)) 
#		flist <- checkforcings(forcings, times, dllname, initforc, 
#			TRUE, fcontrol)
    if (!is.character(func)) 
        stop("`func' must be a *name* referring to a function in a dll")
    if (is.loaded(func, PACKAGE = dllname)) {
        Func <- getNativeSymbolInfo(func, PACKAGE = dllname)$address
    }
    else stop(paste("cannot run DLLfunc: dyn function not loaded: ", 
                        func))
    dy <- rep(0, n)
    storage.mode(y) <- storage.mode(dy) <- "double"
    out <- .Call("call_DLL", y, dy, as.double(times[1]), Func, 
            ModelInit, parms, as.integer(nout), as.double(rpar), 
            as.integer(ipar), as.integer(1), flist, PACKAGE = "deSolve")
    vout <- if (nout > 0) 
                out[(n + 1):(n + nout)]
            else NA
    out <- list(dy = out[1:n], var = vout)
    if (!is.null(Ynames)) 
        names(out$dy) <- Ynames
    if (!is.null(outnames)) 
        names(out$var) <- outnames
    return(out)
}
