#require(testthat)
context("correlations")
# Furher context: fCheckColNames, fCheckColNumeric, fCheckOutsideRange
#library(mvrtnorm)

nParm <- 4L
set.seed(0815)
Sigma <- outer(1:nParm, 1:nParm, .corG, psi=2)
cholSigma <- chol(Sigma)

thetaTrue <- structure( as.vector(rmvnorm(1L, mean=numeric(nParm), sigma=Sigma)), names=paste("a",1:nParm,sep=""))

mod1 <- function(X,mu,Sigma,...){
    dmvnorm(t(X), mean=mu, sigma=Sigma, log=TRUE) 
}
mod1b <- function(
    ### Unscaled Probability density for each column of X of a multivariante distribution 
    X               ##<< numeric matrix: each column a realization of the distribtution
    , mu            ##<< numeric vector: mean of the distribtuion 
    , Sigma         ##<< Covariance matrix
    , cholSigma=chol(Sigma)     ##<< upper triangular Cholesky decomposition of the Covariance matrix, passed for efficiency
){ 
    B <- X-mu
    #Bs <- Sigma^-1 %*% B
    Bs <- backsolve(cholSigma, forwardsolve(t(cholSigma),B))
    sapply( 1:ncol(B), function(i){ B[,i] %*% Bs[,i] })
}

nObs <- 120L
X <- t(rmvnorm(nObs, mean=thetaTrue, sigma=Sigma))
tmp1 <- mod1( X[,1:5], thetaTrue, cholSigma=cholSigma, Sigma=Sigma)
tmp2 <- -0.5*mod1b( X[,1:5], thetaTrue, cholSigma=cholSigma, Sigma=Sigma)
plot(tmp1 ~ tmp2)

#dmvnorm really providing the same information with some additive constant

obsTrue <- mod1(X, thetaTrue, cholSigma=cholSigma, Sigma=Sigma)

plot(density(obsTrue))
LScaled <- -(obsTrue - max(obsTrue))
plot( X[1,], X[2,], col=heat.colors(21)[21-round(LScaled/max(LScaled)*20)])
points( thetaTrue[1], thetaTrue[2], pch="x")


sdObs <- 0.2  #sd(obsTrue)/5
obs <- obsTrue + rnorm( length(obsTrue), sd=sdObs)
#plot( obs ~ obsTrue )

argsMod <- list(
        X=X
        ,Sigma=Sigma
        ,cholSigma=cholSigma
)
argsFLogDensity <- c(argsMod, list(
                obs=obs
                ,sdObs=sdObs
        ))
denNorm <- function(theta, intermediates, logDensityComponents
        , obs, sdObs
        , X, cholSigma
        , ...   ##<< further arguments to mod1, e.g. bias
){
    pred <- mod1(X, mu=theta, cholSigma=cholSigma, ...)
    misfit <- pred - obs
    Sk <- sum( (misfit/sdObs)^2 )
    structure( -1/2 * Sk, names="obs" )
}
logDensityComponents0 <- do.call( denNorm, c( list( thetaTrue, list(), numeric(0) ), argsFLogDensity) )
maxLogDensity <- -1/2*nObs


theta0 <- thetaTrue; theta0[] <- 0
covarTheta0 <- diag(1.5, nrow=length(theta0))

test_that("blockSampling", {
            thetaBlockSpec <- blockSpec(names(theta0),,
                    new("MetropolisBlockUpdater",
                            fLogDensity = denNorm,
                            argsFLogDensity = argsFLogDensity,      
                            logDensityComponentNames = names(logDensityComponents0)
                    ))
            blockSpecs = list( theta=thetaBlockSpec )
            #sampler <- newPopulationSampler( blockSpecs, theta0, covarTheta0, nPopulation=2L )
            #sampler <- setupAndSample(sampler, nSample=120L, thin=4L)
            sampler <- newPopulationSampler( blockSpecs, theta0, covarTheta0, nPopulation=2L )
            sampler <- setupAndSample(sampler, nSample=120L, thin=4L)
            sampleLogs <- getSampleLogs(sampler)
            logsTail <- subsetTail(sampleLogs, 0.4)
            sample1 <-t(getParametersForAllChains(logsTail)) 
            .tmp.f <- function(){
                plot( asMcmc(sampleLogs), smooth=FALSE )
                summary(sample1)
                plot( asMcmc(logsTail), smooth=FALSE )
                .mean
                thetaTrue
                pairs(sample1)
            }
            .mean <- colMeans(sample1)
            .sd <- apply(sample1, 2, sd)
            expect_isOfMagnitude( .sd, 0.02 ) # for N=3
            .pthetaTrue <- pnorm(thetaTrue, mean=.mean, sd=.sd)
            expect_isInInterval( .pthetaTrue )
        })

tmp.inspect.jumps <- function(){
    #recover in PopulationSampler .proposeSampleAndLogRange
    #str(jumps)
    jump <- jumps$jump
    dim(jump) <- c(nrow(jump), prod(dim(jump)[2:3]))
    x <- t(jump)
    pairs(x)
    #ok: indeed the correlations between jumps map the correlations in the sample
    # if a sample is far outside the valley, it will be hard to get back:
    #  correlation structure not good preserved, jumps are small
    # with increasing number of dimensions it gets likeliy that 
    # a sample is outside the valley for one of the multitude of correlations
}
