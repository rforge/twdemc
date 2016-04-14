.tmp.f <- function(){
    library(profr)
    library(plyr)
    library(microbenchmark)
}

#---------------- setup
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


modelPred <- intermediateSpec(
        function(theta, intermediates
                ,xSparse, xRich, ...
        ){
            pred <- modTwTwoDenEx1(theta, xSparse=xSparse, xRich=xRich, ...) 
            return(pred)
        }
        ,argsUpdateFunction=list(xSparse=xSparse, xRich=xRich)
        ,parameters = c("a","b")
)
intermediateSpecs <- list( modelPred = modelPred )
modelPred0 <- getUpdateFunction(modelPred)(theta0, list(), xSparse, xRich) 


denSparse <- function( theta, intermediates, logDensityComponents
        ,obs, sdObs, ...
){
    pred <- intermediates[[1]] 
    misfit <- obs$y1 - pred$y1
    structure( -1/2 * sum((misfit/sdObs$y1)^2), names='y1' )
}
logDensityComponentsSparse0 <- denSparse( theta0, list(modelPred=modelPred0), numeric(0), obs, sdObsTrue )

denRich <- function( theta, intermediates, logDensityComponents
        , obs, sdObs, ...
){
    pred <- intermediates[[1]] 
    misfit <- obs$y2 - pred$y2
    structure( -1/2 * sum((misfit/sdObs$y2)^2), names='y2' )
}
logDensityComponentsRich0 <- denRich( theta0, list(modelPred=modelPred0), numeric(0), obs, sdObsTrue )


argsFLogDensity = list(obs=obs, sdObs=sdObs)

sparseBlockSpec <- blockSpec("a",names(theta0),
        new("MetropolisBlockUpdater",
                fLogDensity = denSparse,
                argsFLogDensity = argsFLogDensity,      
                logDensityComponentNames = names(logDensityComponentsSparse0)
        )
        ,intermediatesUsed=c("modelPred")
)
richBlockSpec <- blockSpec("b",names(theta0),
        new("MetropolisBlockUpdater",
                fLogDensity = denRich,
                argsFLogDensity = argsFLogDensity,      
                logDensityComponentNames = names(logDensityComponentsRich0)
        )
        ,intermediatesUsed=c("modelPred")
)
blockSpecs = list( sparse=sparseBlockSpec, rich=richBlockSpec )

#------------------ profiling chain sampler
blockUpdaters <- newBlockUpdaters(blockSpecs, intermediateSpecs, names(theta0))
thin <- 2L
chainSampler0 <- new("ChainSamplerImpl", blockUpdaters=blockUpdaters, thin=thin)
#        
nGen <- 100L
#nGen <- 500L
steps <- t(rmvnorm( nGen, mean=c(a=0,b=0), sigma=covarTheta0 ))
stepInfoRange <- new("StepInfoRangeImpl", step=steps, nLogDensityComponents=2L)
# starting conditions not in confidence interval 
chainState <- new("ChainState", parameters=c(a=80,b=-20), logDensityComponents=c(y1=NA_real_,y2=NA_real_), intermediates=list(modelPred=list()))
#chainState <- new("ChainState", parameters=fxLinReg$theta0, logDensityComponents=c(obs1=NA_real_))
chainSampler <- setRangeSpecs(chainSampler0,
        chainState=chainState, nInterval=nGen%/%thin, stepInfoRange=stepInfoRange)


.tmp.profileS4Dispatch <- function(){
    # 10 times as fast -> replace most called S4 methods by direct calls
    updater <- getBlockUpdaterForIndex(chainSampler@blockUpdaters, 1L)
    microbenchmark( getBlockIndicesIParametersToUpdate(updater),  getBlockIndicesIParametersToUpdate(updater) )
}

res <- profr(
        sampleRange(chainSampler)
        ,interval = 0.005)
plot(res)



arrange(subset(res, level==10 & start != 0), time)
# long record sample
plot( subset(res, level >= 10 & g_id==10))
plot( subset(res, level >= 15 & g_id==16))
# no  

#------------------ profiling Population sampler
sampler <- newPopulationSampler( blockSpecs, theta0, covarTheta0, nPopulation=1L
        , intermediateSpecifications=intermediateSpecs )

res <- profr( sampler <- setupAndSample(sampler, nSample=200L) )
plot(res)

arrange(subset(res, level==13 & start != 0), time)

# sum execution times 
sumTimes <- ddply( res, .(level, f), function(ds){ c(time=sum(ds$time)) })
head(arrange(sumTimes, rev(time)),20)







