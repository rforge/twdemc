#require(testthat)
context("TwoDenEx1")
# Furher context: fCheckColNames, fCheckColNumeric, fCheckOutsideRange


modTwTwoDenEx1 <- function(
        theta           ##<< model parameters a and b
        , xSparse       ##<< numeric vector of Sparse input/output relationship
        , xRich         ##<< numeric vector of rich input/output relationship
        , thresholdCovar=0      ##<< model structural deficiency
){
    list(  y1 = as.numeric(theta[1]*xSparse + theta[2]*mean(xRich)/10)
           ,y2 = as.numeric(theta[1]*xSparse[1] + theta[2]*pmax(0,xRich-thresholdCovar) ))
}

set.seed(0815)      # for reproducable results
xSparse <- c(1,runif(9,min=0.5,max=1.5))  # only 10 observations
xRich <- runif(1000,min=.7,max=1)  	 # 1000 observations	
thetaTrue <- c( a=1, b=2 )
theta0 <- thetaTrue
covarTheta0 <- diag(0.2*thetaTrue)
obsTrue <- modTwTwoDenEx1( thetaTrue, xSparse=xSparse,xRich=xRich
        , thresholdCovar=0.3)
sdObsTrue <- sdObs <- list( 
        y1=mean(obsTrue$y1)*0.06
        ,y2=mean(obsTrue$y2)*0.02 )
obs <- list(
        y1 = obsTrue$y1 + rnorm(length(obsTrue$y1),sd=sdObsTrue$y1) 
        ,y2 = obsTrue$y2 + rnorm(length(obsTrue$y2),sd=sdObsTrue$y2) )


denSparse <- function( theta, intermediates, logDensityComponents
    ,xSparse, xRich, obs, sdObs, ...
){
    pred <- modTwTwoDenEx1(theta, xSparse=xSparse, xRich=xRich, ...) 
    misfit <- obs$y1 - pred$y1
    structure( -1/2 * sum((misfit/sdObs$y1)^2), names='y1' )
}
logDensityComponentsSparse0 <- denSparse( theta0, list(), numeric(0), xSparse, xRich, obs, sdObsTrue )

denRich <- function( theta, intermediates, logDensityComponents
        ,xSparse, xRich, obs, sdObs, ...
){
    pred <- modTwTwoDenEx1(theta, xSparse=xSparse, xRich=xRich, ...) 
    misfit <- obs$y2 - pred$y2
    structure( -1/2 * sum((misfit/sdObs$y2)^2), names='y2' )
}
logDensityComponentsRich0 <- denRich( theta0, list(), numeric(0), xSparse, xRich, obs, sdObsTrue )


argsFLogDensity = list(xSparse=xSparse, xRich=xRich, obs=obs, sdObs=sdObs)

# with original unbiased thresholdCovar
#argsFLogDensity = list(xSparse=xSparse, xRich=xRich, obs=obs, sdObs=sdObs, thresholdCovar = 0.3)


test_that("blockSampling", {
            sparseBlockSpec <- blockSpec("a",names(theta0),
                    new("MetropolisBlockUpdater",
                            fLogDensity = denSparse,
                            argsFLogDensity = argsFLogDensity,      
                            logDensityComponentNames = names(logDensityComponentsSparse0)
                    )
            )
            richBlockSpec <- blockSpec("b",names(theta0),
                    new("MetropolisBlockUpdater",
                            fLogDensity = denRich,
                            argsFLogDensity = argsFLogDensity,      
                            logDensityComponentNames = names(logDensityComponentsRich0)
                    )
            )
            blockSpecs = list( sparse=sparseBlockSpec, rich=richBlockSpec )
            
            sampler <- newPopulationSampler( blockSpecs, theta0, covarTheta0, nPopulation=1L )
            sampler <- setupAndSample(sampler, nSample=60L)
            sampleLogs <- getSampleLogs(sampler)
            plot( asMcmc(sampleLogs), smooth=FALSE )
            
            stacked <- stackChainsInPopulation(subsetTail(getSampleLogs(sampler),0.8))
            sample1 <- t(getParametersForPopulation(stacked, 1L)[,,1L])
            summary(sample1)
            
            plot( asMcmc(stacked), smooth=FALSE )
            .mean <- colMeans(sample1)
            .sd <- apply(sample1, 2, sd)
            expect_isOfMagnitude( c(a=0.04,b=0.04), .sd )
            .pthetaTrue <- pnorm(c(a=1.04,b=1.2), mean=.mean, sd=.sd)
            expect_isInInterval( .pthetaTrue ) 
            
            ss <- cbind( sample1,  t(getLogDensitiesForPopulation(stacked,1L)[,,1L]))
            r <- computeSampleRanksForPopulation(stacked, 1L)[,1L]
        })







    