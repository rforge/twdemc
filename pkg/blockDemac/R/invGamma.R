#' @export
getInvGammaParsFromMeanAndMode <- function(
        ### calculate parameters of inverse gamma function based on given mean and mode
        mean    ##<< mean of the distribution
        ,mode   ##<< mode of the distribution
){
    # https://en.wikipedia.org/wiki/Inverse-gamma_distribution
    alpha <- as.numeric((mean+mode)/(mean-mode))
    beta <- as.numeric(mean*(alpha-1))
    ##value<< numeric vector with components
    c( alpha=alpha      ##<< shape parameter
            , beta=beta )    ##<< scale parameter
}

#' @export
rInvGamma <- function(
    ### Density function of the inverse gamma distribution. 
    n           ##<< Number of draws from the distribution.
    , shape     ##<< Scalar shape parameter.
    , scale = 1 ##<< Scalar scale parameter (default value one).
){
    # originally from the MCMCpack package. But this has too many dependencies
    # e.g. graph from Bioconductor and installation problems
    #
    ##details<<
    ## An inverse gamma random variable with shape a and scale b 
    ## has mean b/(a-1) (assuming a>1) and variance (b^2)/((a-1)^2 (a-2)) (assuming a>2).
    ## The parameterization is consistent with the Gamma Distribution in the stats package
    #
    ##references<<
    ## Andrew Gelman, John B. Carlin, Hal S. Stern, and Donald B. Rubin. 2004. Bayesian Data Analysis. 2nd Edition. Boca Raton: Chapman & Hall.
    #
    ##seealso<<
    ## \code{\link{GammaDist}}
    ##value<<
    ## n draws from the inverse Gamma distribution. 
    return(1/rgamma(n = n, shape = shape, rate = scale))
}
attr(rInvGamma,"ex") <- function(){
    draws <- rInvGamma(10, 3.2)    
}

##' @export
#rMvNorm <- function(
#        ### call to \code{\link{rmvnorm}} from package mvtnorm
#        ...
#){ 
#    ##details<< 
#    ## Functions from other packages are not visible in the workspace
#    ## However, need to call it in a function put to the cluster, without requireing
#    ## the user to explicitely load package on the cluster.
#    rmvnorm(...) 
#}

#' @export
getGammaParsFromMeanAndVariance <- function(
        ### calculate parameters of  gamma function based on given mean and variance
        mean    ##<< mean of the distribution
        ,var   ##<< variance of the distribution
){
    # https://en.wikipedia.org/wiki/Gamma_distribution
    theta <- var/mean
    k <- mean/theta
    ##value<< numeric vector with components
    c( k=k      ##<< shape parameter
      ,theta=theta )    ##<< scale parameter
}
attr(getGammaParsFromMeanAndVariance,"ex") <- function(){
    m <- 1
    v <- 0.2
    #m <- 1/3; v=(1/3.2)^2
    parsGamma <- getGammaParsFromMeanAndVariance(m,v)
    parsGamma["k"] * parsGamma["theta"] - m < 1e10
    parsGamma["k"] * parsGamma["theta"]^2 -v < 1e10
    tmp <- exp(seq(-10,log(m+2*sqrt(v)),length.out=91))
    plot( dgamma(tmp, parsGamma["k"], scale=parsGamma["theta"]) ~ tmp); abline(v=m)
}




