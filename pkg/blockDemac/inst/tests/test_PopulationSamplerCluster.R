#require(testthat)
context("PopulationSamplerCluster")


set.seed(0815)      # for reproducible results
simpleModel <- function(theta,xval){ theta[1]*xval }

sdObsTrue <- 0.01
thetaTrue <- c(theta=0.65) # the true parameter vector

nObs <- 11
x <- seq(0,4,length.out=nObs+1)[-1]
obsTrue <- simpleModel(thetaTrue,x)/(1+x/20)
obs <- obsTrue + rnorm(length(x), sd=sdObsTrue)
sd2Obs <- sdObsTrue^2

thetaHat <- structure(coef(lm(obs ~x-1)), names="theta")

theta0 = c(theta=as.numeric(thetaHat), logPsi=log(0.6), logSd2Discr=log(0.5))
covarTheta0 = diag(c(theta0*0.2)^2, nrow=length(theta0))    # spread for initial population 
streamSuffix <- "obs1"
theta0N <- theta0; names(theta0N)[2:3] <- paste(names(theta0)[2:3],streamSuffix,sep="_")  



predProcSpec <- intermediateSpec(
        #dumpOnError(
            function(theta, intermediates, ..., xVal, fSimpleModel ){
                fSimpleModel(theta["theta"],xVal)
            }
        #)
        # additional arguments and functions called inside update function 
        # will not be available in remote process.
        # need to either export before or pass them as arguments and list in argsUpdateFunction
        # different names are required x=x will not work (why?)
        ,argsUpdateFunction=list(xVal=x, fSimpleModel=simpleModel)        
        ,parameters = c("theta")
        ,intermediates = character(0)  
)

# Metropolis block sampling theta, based on intermediate predProcSpec
logDenTheta <- function( theta, intermediates, logDensityComponents, obs, sd2Obs ){
    predProc <- intermediates$predProc_obs1 # heed the streamSuffix
    residuals <- obs - (predProc)
    return( c(obs1= -1/2 * sum( ((residuals)^2)/sd2Obs )) )
}
testLogDenTheta <- logDenTheta(theta0, list(predProc_obs1=simpleModel(thetaHat,x)), numeric(1), obs=obs, sd2Obs=sd2Obs )
logDenThetaSpec <- blockSpec("theta","theta",    
        new("MetropolisBlockUpdater",
                fLogDensity=logDenTheta,
                argsFLogDensity=list( obs=obs, sd2Obs=sd2Obs),
                logDensityComponentNames = names(testLogDenTheta)
        )
        ,intermediatesUsed=paste(c("predProc"),streamSuffix,sep="_")
)
#
test_that("sampling in non-parallel setting",{
            sampler <- newPopulationSampler( 
                    list(bTheta=logDenThetaSpec)
                    , theta0N[1], covarTheta0[1,1, drop=FALSE]
                    ,intermediateSpecifications=list(predProc_obs1=predProcSpec)
            )
            
            sampler <- suppressWarnings(setupAndSample(sampler, nSample=3L))
    sampleLogs <- getSampleLogs(sampler)
    #plot( asMcmc(sampleLogs), smooth=FALSE )
})

test_that("sampling in parallel setting",{
    # note that this uses the installed package, loaded to the cluster
    # maybe need to do genDep() and instPkg() before this test works
    require(parallel)
    try( stopCluster(cl), silent=TRUE )
    cl <- makeCluster(2L)
    on.exit(stopCluster(cl))
    #clusterExport(cl, c("mcCommon", ".sampleRangeOneChain"), envir=environment())
    clusterEvalQ(cl, library(blockDemac) )
    #clusterEvalQ(cl, library(twMisc) )     # dumpOnError
    clusterExport(cl, c("simpleModel"), envir=environment() ) #need envir argument, because test runs not in global environment
    
    sampler <- newPopulationSampler( 
            list(bTheta=logDenThetaSpec)
            , theta0N[1], covarTheta0[1,1, drop=FALSE]
            ,intermediateSpecifications=list(predProc_obs1=predProcSpec)
    )
    sampler@cluster <- cl
    sampler <- suppressWarnings(setupAndSample(sampler, nSample=3L))
    sampleLogs <- getSampleLogs(sampler)
    #plot( asMcmc(sampleLogs), smooth=FALSE )
})


