linearModelFunctionExample <- function(
    ### example model function: y=a+bx
    theta,	##<< parameter vector with names a and b
    xval	##<< additional argument aside parameter vector: numeric vector of a single covariate
){ 
    theta["a"] + theta["b"]*xval 
}


