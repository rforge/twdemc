#require(testthat)
context("denNormal15")
# Furher context: fCheckColNames, fCheckColNumeric, fCheckOutsideRange

nParm <- 15L
set.seed(0815)
thetaTrue <- structure( rnorm(nParm), names=paste("a",1:nParm,sep=""))
theta0 <- thetaTrue; theta0["a1"] <- 0
covarTheta0 <- diag(0.5, nrow=length(theta0))
mod1 <- function(theta){ (theta) }

obsTrue <- mod1(thetaTrue)
sdObs <- 0.1
obs <- obsTrue + rnorm( length(obsTrue), sd=sdObs)


denNorm <- function(theta, intermediates, logDensityComponents, obs=obs, sdObs=sdObs){
    pred <- theta   #mod1(theta) 
    misfit <- pred - obs
    Sk <- sum(misfit^2)
    structure( -1/2 * sum((misfit/sdObs)^2), names='obs' )
}
logDensityComponents0 <- denNorm( theta0, list(), numeric(0), obs=obs, sdObs=sdObs )

argsFLogDensity = list(obs=obs, sdObs=sdObs)

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
            sampler <- newBatchSampler( blockSpecs, theta0, covarTheta0, nPopulation=2L, thin=4L )
            sampler <- suppressMessages(
                    #p1 <- profr(
                    sampler <- sampleBatches(sampler, nBatchMax=12L )
                    #)
            )
            sampleLogs <- getSampleLogs(sampler)
            logsTail <- subsetTail(sampleLogs, 0.5)
            sample1 <-t(getParametersForAllChains(sampleLogs)) 
            .tmp.f <- function(){
                getBatchLog(sampler)
                plot( asMcmc(sampleLogs), smooth=FALSE )
                summary(sample1)
                logsTail <- subsetTail(sampleLogs, nSample=120)
                plot( asMcmc(logsTail), smooth=FALSE )
            }
            .mean <- colMeans(sample1)
            .sd <- apply(sample1, 2, sd)
            expect_isOfMagnitude( .sd, sdObs )
            .pthetaTrue <- pnorm(thetaTrue, mean=.mean, sd=.sd)
            expect_isInInterval( .pthetaTrue ) 
        })

.tmp.f <- function(){
    sampler <- newBatchSampler( blockSpecs, theta0, covarTheta0, nPopulation=2L, thin=4L )
    sampler <- sampleBatches(sampler, nBatchMax=12L )
    
}

.testParallel <- function(){
    require(parallel)
    try( stopCluster(cl), silent=TRUE )
    cl <- makeCluster(2L)
    on.exit(stopCluster(cl))
#clusterExport(cl, c("mcCommon", ".sampleRangeOneChain"), envir=environment())
    clusterEvalQ(cl, library(blockDemac) )
    #clusterEvalQ(cl, library(twMisc) )     # dumpOnError
    
    sampler@cluster <- cl
}








    