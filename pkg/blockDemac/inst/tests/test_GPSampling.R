#require(testthat)
context("GPSampling")

# see inst/TowDenGP/TwoDenGP.R
modTwTwoDenEx1 <- function(
        theta           ##<< model parameters a and b
        , xSparse       ##<< numeric vector of Sparse input/output relationship
        , xRich         ##<< numeric vector of rich input/output relationship
        , thresholdCovar=0      ##<< model structural deficiency
){
    list(  y1 = as.numeric(theta[1]*xSparse + theta[2]*mean(xRich)/10)
            ,y2 = as.numeric(theta[1]*xSparse[1] + theta[2]*pmax(0,xRich-thresholdCovar) ))
}

xSparse <- sort(c(1,runif(9,min=0.5,max=1.5)))  # only 10 observations
xRich <- sort(runif(1000,min=.7,max=1))    # 1000 observations    
# for plotting predictions, define a grid
xGridRich <- seq( xRich[1], xRich[length(xRich)], length.out=30)
xGridSparse <- seq( xSparse[1], xSparse[length(xSparse)], length.out=30)

thetaTrue <- c( a=1, b=2 )
theta0 <- thetaTrue
covarTheta0 <- diag(0.2*thetaTrue)
thresholdCovarTrue=0.3
obsTrue <- modTwTwoDenEx1( thetaTrue, xSparse=xSparse,xRich=xRich
        , thresholdCovar=thresholdCovarTrue)
sdObsTrue <- sdObs <- list( 
        #y1=mean(obsTrue$y1)*0.06
        #,y2=mean(obsTrue$y2)*0.02 )
        #y1=mean(obsTrue$y1)*0.04
        #,y2=mean(obsTrue$y2)*0.03 )
        y1=mean(obsTrue$y1)*0.05
        ,y2=mean(obsTrue$y2)*0.03 )
obs <- list(
        y1 = obsTrue$y1 + rnorm(length(obsTrue$y1),sd=sdObsTrue$y1) 
        ,y2 = obsTrue$y2 + rnorm(length(obsTrue$y2),sd=sdObsTrue$y2) )

# locations where model+discrepancy will be sampled        
xPredSparse <- c(range(xSparse), 1.0)   
xPredRich <- c(range(xRich), 0.8)

thresholdCovarBiased <- 0.1


predSpec <- intermediateSpec(
        function(theta, intermediates, xSparseVal, xRichVal, ... ){
            modTwTwoDenEx1(theta, xSparse=xSparseVal, xRich=xRichVal, ...) 
        }
        ,argsUpdateFunction=list(xSparseVal=c(xSparse,xPredSparse), xRichVal=c(xRich,xPredRich)
                , thresholdCovar=thresholdCovarBiased)
        ,parameters = c("a","b")
        ,intermediates = character(0)  
)
pred0 <- getUpdateFunction(predSpec)( theta0, list(), xSparse, xRich )


xSparseRange <- diff(range(xSparse))
blocksSparse <- compileDiscrepancyBlocksForStream(x=xSparse, obs=obs$y1, sd2Obs=(sdObs$y1)^2
        , predProcSpec = list(pred=predSpec), predProcSpecEntry="y1"
        , streamSuffix="sparse", xPred=xPredSparse
        ,priorGammaPsiPars = getGammaParsFromMeanAndVariance(xSparseRange/3, (xSparseRange/3.2)^2 )
)

xRichRange <- diff(range(xRich))
blocksRich <- compileDiscrepancyBlocksForStream(x=xRich, obs=obs$y2, sd2Obs=(sdObs$y2)^2
        , predProcSpec = list(pred=predSpec), predProcSpecEntry="y2"
        , streamSuffix="rich", xPred=xPredRich
        ,priorGammaPsiPars = getGammaParsFromMeanAndVariance(xRichRange/3, (xRichRange/3.2)^2)
)

namesPred <- c(blocksSparse$namesProcessSamples, blocksRich$namesProcessSamples )
thetaPred <- structure( rep(1,length(namesPred)), names=namesPred)
theta0Md <- c(theta0, thetaPred
        , logPsi_sparse=0.2, logSd2Discr_sparse=0.5
        , logPsi_rich=0.2, logSd2Discr_rich=0.5)
covarTheta0Md <- diag((0.2*theta0Md)^2)


logDenThetaMd <- function( theta, intermediates, logDensityComponents
        ,xSparse, xRich, obs, sdObs
        , priorDiscrFacSparse= length(xSparse)/2
        , priorDiscrFacRich= length(xRich)/2
        , ...
){
    intSparse <- intermediates[["deltaA_sparse"]]
    intRich <- intermediates[["deltaA_rich"]]
    deltaSparse <- intSparse$deltaA
    deltaRich <- intRich$deltaA
    #misfit1 <- obs$y1 - (intermediates[["pred"]]$y1[seq_along(xSparse)] + deltaSparse)
    #misfit2 <- obs$y2 - (intermediates[["pred"]]$y2[seq_along(xRich)] + deltaRich)
    #logPriorSd2DiscrSparse <- priorDiscrFacSparse*intermediates$deltaA_sparse$logPriorInvGammaSd2Discr 
    #logPriorSd2DiscrRich <- priorDiscrFacRich*intermediates$deltaA_rich$logPriorInvGammaSd2Discr
    logPriorSd2DiscrSparse <- intermediates$deltaA_sparse$logPDeltaO  
    logPriorSd2DiscrRich <- intermediates$deltaA_rich$logPDeltaO
    #logPObs <- structure(-1/2*c(sum((misfit1/sdObs$y1)^2), sum((misfit2/sdObs$y2)^2)), names=c('ySparse','yRich') )   
    logPObs <- structure(-1/2*c(intSparse$SSQpred, intRich$SSQpred), names=c('ySparse','yRich') )   
    #if( logPObs["yRich"] < -515 ) recover()
    res <- c( logPObs, mdSparse=logPriorSd2DiscrSparse, mdRich=logPriorSd2DiscrRich)
    res
}

resLogDenThetaMd0 <- logDenThetaMd( theta0, list(pred=pred0
                ,deltaA_rich=list(deltaA=0, logPriorInvGammaSd2Discr=-2, logPDeltaO=-2, SSQpred=13.9)
                , deltaA_sparse=list(deltaA=0, logPriorInvGammaSd2Discr=-2, logPDeltaO=-2, SSQpred=13.9))
        , xSparse=xSparse, xRich=xRich, obs=obs, sdObs=sdObs )

bThetaMdSpec <- blockSpec(
        c("a","b"),,    
        new("MetropolisBlockUpdater",
                #logDenThetaMd,
                fLogDensity=function(...){ logDenThetaMd(...) },    # allows tracing
                argsFLogDensity=list( xSparse=xSparse, xRich=xRich, obs=obs, sdObs=sdObs),
                logDensityComponentNames = names(resLogDenThetaMd0)
        )
        ,intermediatesUsed=c("pred", "deltaA_sparse", "deltaA_rich")
)

cSubSpace <- blocksRich$fConstrainingSubSpace( 
        blocksSparse$fConstrainingSubSpace(new("SubSpace")) )
# know that correlation length is in the magnitude of the locations    
cSubSpace@lowerParBounds["logPsi_rich"] <- log(diff(range(xRich))/5)      


test_that("sampling in non-parallel setting",{
            sampler <- newPopulationSampler( 
                    blockSpecifications=c(list(bTheta=bThetaMdSpec)
                            ,blocksSparse$blockSpecs, blocksRich$blockSpecs)
                    , theta0Md, covarTheta0Md
                    , nPopulation=1L
                    ,intermediateSpecifications=c(blocksSparse$intermediateSpecs
                            , blocksRich$intermediateSpecs[-1L])   #-1L to not include the prediction intermediate twice 
                    ,subSpace=cSubSpace
            )
            sampler <- setupAndSample(sampler, nSample=4L)
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
    clusterExport(cl, c("modTwTwoDenEx1","logDenThetaMd"), envir=environment() )
    
    sampler <- newPopulationSampler( 
            blockSpecifications=c(list(bTheta=bThetaMdSpec)
                    ,blocksSparse$blockSpecs, blocksRich$blockSpecs)
            , theta0Md, covarTheta0Md
            , nPopulation=1L
            ,intermediateSpecifications=c(blocksSparse$intermediateSpecs
                    , blocksRich$intermediateSpecs[-1L])   #-1L to not include the prediction intermediate twice 
            ,subSpace=cSubSpace
    )
    sampler@cluster <- cl
    sampler <- setupAndSample(sampler, nSample=4L)
    sampleLogs <- getSampleLogs(sampler)
    #plot( asMcmc(sampleLogs), smooth=FALSE )
})


