#require(testthat)
#context("denNormal15")
# Furher context: fCheckColNames, fCheckColNumeric, fCheckOutsideRange
#library(mvtnorm)

#' @export
.compileMultivariateNormalProblem <- function(
        nParm=15L
        ,Sigma=diag(nrow=nParm)     ##<< Covariance matrix of the parameters
        ,nObs = c(rich=600L, sparse=80L) 
        , bias = c(rich=-1)         ##<< numeric vector: For each stream a Bias scalar that is added to the middle component of theta
            ##<< if it is a scalar, then it applies to the first stream and the other streams have zero bias
){
    set.seed(0815)
    thetaTrue <- structure( as.vector(rmvnorm(1L, rep(0,nParm), sigma=Sigma)), names=paste("a",1:nParm,sep=""))
    cholSigma <- chol(Sigma)
    #
    mod1 <- function(
            theta=thetaTrue ##<< the parameters
            , Xlist         ##<< list with entry for each data stream: matrix with rows of covariates for which to calculate multivariate normal
            , cholSigma=chol(Sigma)    ##<< Cholesky decomposition of covariance matrix, usually passed for efficiency
            , Sigma         ## covariance matrix
            , bias=rep(0,length(Xlist)) ##<< numeric vector: For each stream a Bias scalar that is added to first half of components of theta
    ){
        thetaList <- rep( list(theta), length(Xlist) )
        #iCompsBiased <- 1:as.integer(length(thetaTrue)/4)   # cannot find suitable parameters
        #iCompsBiased <- seq_along(thetaTrue)[-1]
        iCompsBiased <- TRUE    # all components
        structure(lapply(seq_along(Xlist), function(iStream){
                            thetaStream <- thetaList[[iStream]]
                            thetaStream[iCompsBiased] <- thetaStream[iCompsBiased] + bias[iStream]
                            #dmvnorm(Xlist[[i]],thetaStream, Sigma, log=TRUE)
                            dmvnormOpt(Xlist[[iStream]],thetaStream, cholSigma=cholSigma, log=TRUE)
                        }),names=names(Xlist))
    }
    #
    Xlist <- lapply( nObs, rmvnorm, mean=thetaTrue, sigma=Sigma)
    obsTrue <- mod1(thetaTrue,Xlist, cholSigma=cholSigma, bias=rep(0, length(nObs)) )
    sdObs <- lapply(obsTrue, function(obsStreamTrue){
                sd(obsTrue[[1]])/5
            })
    obs <- structure(lapply(names(obsTrue), function(stream){
                        obsTrue[[stream]] + rnorm( length(obsTrue[[stream]]), sd=sdObs[[stream]])
                    }), names=names(obsTrue))
#plot( obs[[1]] ~ obsTrue[[1]] )
    #
    biasFormatted = if( length(bias)==1L ){
       biasFormatted <- numeric(length(nObs)); biasFormatted[1] <- bias; biasFormatted
    } else if( length(bias) != length(nObs) ) {
        stop("components and order in bias must match those of nObs.")
    } else bias
    #    
    argsMod <- list(
            Xlist=Xlist
            #,Sigma=Sigma
            ,cholSigma=cholSigma
            #,bias=-2
            ,bias=biasFormatted
    )
    argsFLogDensity <- c(argsMod, list(
                    obs=obs
                    ,sdObs=sdObs
            ))
    #
    denNorm <- function(theta, intermediates, logDensityComponents
            , obs, sdObs
            , Xlist, cholSigma
            , ...   ##<< further arguments to mod1, e.g. bias
    ){
        pred <- mod1(theta, Xlist,  cholSigma=cholSigma, ...) 
        SkStreams <- sapply( names(pred), function(stream){
                    misfit <- pred[[stream]] - obs[[stream]]
                    Sk <- sum( (misfit/sdObs[[stream]])^2 )
                })
        structure( -1/2 * SkStreams, names=names(pred) )
    }
    logDensityComponents0 <- do.call( denNorm, c( list( thetaTrue, list(), numeric(0) ), argsFLogDensity) )
    maxLogDensity <- -1/2*nObs
    #
    thetaBlockSpec <- blockSpec(names(thetaTrue),,
            new("MetropolisBlockUpdater",
                    fLogDensity = denNorm,
                    argsFLogDensity = argsFLogDensity,      
                    logDensityComponentNames = names(logDensityComponents0)
                    ,maxLogDensity = maxLogDensity
            ))
    #        
    ans <- list(
            nParm=nParm
            ,sigma=Sigma
            ,nObs=nObs
            ,thetaTrue=thetaTrue
            ,covarTheta0=1.5*Sigma
            ,bias=biasFormatted
            ,model=mod1
            ,Xlist=Xlist
            ,obsTrue=obsTrue
            ,obs=obs
            ,sdObs=sdObs
            ,argsMod=argsMod
            ,argsFLogDensity=argsFLogDensity
            ,fLogDen=denNorm
            ,logDensityComponents0=logDensityComponents0
            ,maxLogDensity=maxLogDensity
            ,thetaBlockSpec=thetaBlockSpec
    )
}

#' @export
dmvnormOpt <- function(
        ### optimized version of dmvnorm, with providing cholesky decomposition
        x, mean, sigma, cholSigma=chol(sigma), log = FALSE  
){
    ##details 
    ## omits argument checking: x must be a matrix with each row quantile vector
    p <- ncol(x)
    tmp <- backsolve(cholSigma, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(cholSigma))) - 0.5 * p * log(2 * 
                    pi) - 0.5 * rss
    names(logretval) <- rownames(x)
    if (log) 
        logretval
    else exp(logretval)
}