.tmp.f <- function(){
    library(RBenchmark)
    library(profr)
    library(plyr)
}

fxLinReg <- .fixtureLinReg1()
lmDummy <- with(fxLinReg, lm( obs ~ xval, weights=1/sdObs^2))		# results without priors
(.expTheta <- coef(lmDummy))
(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=names(fxLinReg$theta0)) )
#
blockUpdaters <- newBlockUpdaters(fxLinReg$blocks, list(), names(fxLinReg$theta0))
thin <- 2L
chainSampler0 <- new("ChainSamplerImpl", blockUpdaters=blockUpdaters, thin=thin)
#        
nGen <- 100L
steps <- t(rmvnorm( nGen, mean=c(a=0,b=0), sigma=.expCovTheta*1.2 ))
stepInfoRange <- new("StepInfoRangeImpl", step=steps)
# starting conditions not in confidence interval 
chainState <- new("ChainState", parameters=c(a=80,b=-20), logDensityComponents=c(obs1=NA_real_))
#chainState <- new("ChainState", parameters=fxLinReg$theta0, logDensityComponents=c(obs1=NA_real_))
chainSampler <- setRangeSpecs(chainSampler0,
        chainState=chainState, nInterval=nGen%/%thin, stepInfoRange=stepInfoRange)



res <- profr(
        sampleRange(chainSampler)
,interval = 0.005)
plot(res)

arrange(subset(res, level==10 & start != 0), time)
# long record sample
plot( subset(res, level >= 10 & g_id==3))
plot( subset(res, level >= 10 & g_id==1))


plot(
        #arrange(subset(res, level >= 10 & (g_id %in% c(2))), time) # look into isAccepteBlocks -> blockUpdaterApply
        subset(res, level >= 10 & (g_id %in% c(4)) ) # look into isAccepteBlocks -> blockUpdaterApply
)   

