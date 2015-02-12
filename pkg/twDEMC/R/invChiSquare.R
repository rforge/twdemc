# used in estimating variance of observations, Gelman et al., 2003, p.51
# see inst/SigmaBlockEx1

updateSigmaByGibbsSamplingInvchisq <- function( 
        ### update Variance by sampling from a scaled inverse Chi-square distribution with using prior information on sigma 
        theta				##<< numeric vector: current state that is used in density function of the block
        ,argsFUpdateBlock=NULL   ##<< list with possible component \code{intermediate} that will be provided to \code{fResid}
        ,fResid             ##<< function to calculate residuals between prediction and observations, must return a numeric vector
        ,nu0=0              ##<< weight of prior: equivalent number of observations, defaults to no prior (zero obs) 
        ,sigma20            ##<< prior estimate of variance, needs to be specified only when giving  nu0 different from zero
        ,...                ##<< further arguemnts to fResid
){
    ##details<< 
    ## Gelman et al., 2003, p.51: Estimating variance with known mean.
    ##
    ## See the vignette on unknown observation uncertainty.
    resid <- if( length(argsFUpdateBlock$intermediate) ){
        resid <- fResid(theta, intermediate=argsFUpdateBlock$intermediate, ...)
    }else{
        resid <- fResid(theta,...)
    }
    n <- length(resid)
    v <- sum( resid^2 ) / n
    if( !length(nu0) || nu0 == 0 ){
        sigma2 <- rinvchisq( 1, n, v )
    }else{
        sigma2 <- rinvchisq( 1, nu0 + n, (nu0*sigma20 + n*v)/(nu0+n) )
    }
    ##value<< list with components
    list(	##describe<<
            accepted=1			##<< boolean scalar: if step was accepted
            , xC=log(sigma2)	##<< numeric vector: components of position in parameter space that are being updated
            , intermediate=attr(resid,"intermediate")  ##<< new intermediate state, if provided by fResid
    )	##end<<
}



dinvchisq <- function (
        ### density function and for the (scaled) inverse-chi-squared distribution.
        x					##<< vector of quantiles
        , df				##<< degrees of freedom parameter, usually represented as nu
        , scale = 1/df		##<< scale parameter, usually represented as lambda.
        , log = FALSE		##<< Logical. If log=TRUE, then the logarithm of the density is returned.
){
    # adopted from copyLeft pakcage LaplaceDemon, as not available for R 2.10 which runs on the cluster
    ##seealso<< \code{\link{dinvchisq}}
    x <- as.vector(x)
    df <- as.vector(df)
    scale <- as.vector(scale)
    if (any(x <= 0)) 
        stop("x must be positive.")
    if (any(df <= 0)) 
        stop("The df parameter must be positive.")
    if (any(scale <= 0)) 
        stop("The scale parameter must be positive.")
    NN <- max(length(x), length(df), length(scale))
    x <- rep(x, len = NN)
    df <- rep(df, len = NN)
    scale <- rep(scale, len = NN)
    nu <- df/2
    dens <- nu * log(nu) - log(gamma(nu)) + nu * log(scale) - 
            (nu + 1) * log(x) - (nu * scale/x)
    if (log == FALSE) 
        dens <- exp(dens)
    dens
}
attr(dinvchisq,"ex") <- function(){
    x <- dinvchisq(1,1,1)
    x <- rinvchisq(10,1)
    
    #Plot Probability Functions
    x <- seq(from=0.1, to=5, by=0.01)
    plot(x, dinvchisq(x,0.5,1), ylim=c(0,1), type="l", main="Probability Function",
            ylab="density", col="red")
    lines(x, dinvchisq(x,1,1), type="l", col="green")
    lines(x, dinvchisq(x,5,1), type="l", col="blue")
    legend(3, 0.9, expression(paste(nu==0.5, ", ", lambda==1),
                    paste(nu==1, ", ", lambda==1), paste(nu==5, ", ", lambda==1)),
            lty=c(1,1,1), col=c("red","green","blue"))	
}

rinvchisq <- function (
        ### random number function and for the (scaled) inverse-chi-squared distribution.
        n				## the number of observations. If length(n) > 1, then the length is taken to be the number required.
        , df
        , scale = 1/df
){
    ##seealso<< \code{\link{rinvchisq}}
    df <- rep(df, len = n)
    scale <- rep(scale, len = n)
    if (any(df <= 0)) 
        stop("The df parameter must be positive.")
    if (any(scale <= 0)) 
        stop("The scale parameter must be positive.")
    x <- (df * scale)/rchisq(n, df = df)
    return(x)
}

