context("den2dCor")

.tmp.f <- function(){
    gridx <- a <- seq(-0.5,2,length.out=91)
    gridy <- seq(-20,+40,length.out=91)
    gridX <- expand.grid(gridx, gridy)
    den2dCor(c(0.8,0.8))
    luden <- apply( gridX, 1, den2dCor ) 
    mLuden <- matrix(luden,nrow=length(gridx))
    #plot( rowSums(mLuden) ~ gridx )
    imax <- which( matrix(luden,nrow=length(gridx))==max(luden), arr.ind=TRUE)
    #c( gridx[ imax[1] ], gridy[ imax[2] ] )
    
    image( gridx, gridy,  mLuden, col = rev(heat.colors(100)), xlab="a", ylab="b" )
    xyMax <- c(x=gridx[ imax[1] ], y=gridy[ imax[2] ])
    
    image( gridx, gridy,  matrix(exp(luden),nrow=length(gridx)), col = rev(heat.colors(100)), xlab="a", ylab="b" )
    points( gridx[ imax[1] ], gridy[ imax[2] ]  )
    
    q20 <- quantile(luden,0.2) 
    plot(density(luden[luden>q20]))
    head(sort(luden,dec=TRUE))
}


.tmp.f <- function(){
    ld1 <- den2dCor(c(0.8,0.8))
    muA1 <- -1.5
    sigmaA1 <- 1.2
    blocks <- list(
            m1 = blockSpec(,,   # empty defaults to all parameters 
                    new("MetropolisBlockUpdater",
                            fLogDensity=den2dCor
                            ,argsFLogDensity=list(muA=muA1, sigmaA=sigmaA1)
                            ,logDensityComponentNames = names(ld1)
                    )
            )
    )  
	theta0 <- c(a=0,b=0) 
	covarTheta0 <- diag(c(a=2,b=2)) 		
    sampler0 <- newPopulationSampler( blocks, theta0, covarTheta0, nPopulation=2L, nChainInPopulation=9L, nChainGroupsInPopulation=3L )
    sampler <- setupAndSample(sampler0, nSample=40L)
    
    sampler <- setupAndSample(sampler, nSample=80L)
    
    sampleLogsTail <- subsetTail(getSampleLogs(sampler),0.6)
    #plot( asMcmc(sampleLogsTail), smooth=FALSE )
    aC <- getAcceptanceTracker(sampler)
    computePopulationAcceptanceRates(aC)
    #
    stacked <- stackChainsInPopulation(sampleLogsTail)    
    ss <- ss0 <- t(getParametersForAllChains(sampleLogsTail))    
    
    image( gridx, gridy,  matrix(exp(luden),nrow=length(gridx)), col = rev(heat.colors(100)), xlab="a", ylab="b" )
    points( gridx[ imax[1] ], gridy[ imax[2] ], pch="X"  )
    points( ss[,1], ss[,2] )
    
    plot( ss[,1], ss[,2], ylim=c(-40,80) )
    plot( density(ss[,1]) )
    
    # here with no groupings
    sampler0UnGrouped <- newPopulationSampler( blocks, theta0, covarTheta0, nPopulation=2L, nChainInPopulation=3L )
    samplerUnGrouped <- setupAndSample(sampler0UnGrouped, nSample=40L*3L)
    
    samplerUnGrouped <- setupAndSample(samplerUnGrouped, nSample=80L*3L)
    
    sampleLogsTailUnGrouped <- subsetTail(getSampleLogs(samplerUnGrouped),0.6)
    #plot( asMcmc(sampleLogsTailUnGrouped), smooth=FALSE )
    aCUngrouped <- getAcceptanceTracker(samplerUnGrouped)
    computePopulationAcceptanceRates(aCUngrouped)
    matplot(getAcceptanceLog(aC)[1L,,], type="l")
    matplot(getAcceptanceLog(aCUngrouped)[1L,,], type="l")
    plot(hist(getAcceptanceLog(aC)[1L,,]))
    plot(hist(getAcceptanceLog(aCUngrouped)[1L,,]))
    
    #
    ssUngrouped <- ss0 <- t(getParametersForAllChains(sampleLogsTailUnGrouped))    
    
    image( gridx, gridy,  matrix(exp(luden),nrow=length(gridx)), col = rev(heat.colors(100)), xlab="a", ylab="b" )
    points( gridx[ imax[1] ], gridy[ imax[2] ], pch="X"  )
    points( ssUngrouped[,1], ssUngrouped[,2] )
    
    plot( ssUngrouped[,1], ssUngrouped[,2], ylim=c(-40,80) )
    plot( density(ssUngrouped[,1]) )
    lines( density(ss[,1]), col="blue")
    gridx2 <- a <- seq(-4,2,,length.out=91)
    pX <- dnorm(gridx2, muA, sigmaA)
    lines( pX ~ gridx2, col="orange")
}

.tmp.denBanana <- function(){
    mleA <- 0.5
    coefA <- twCoefLogitnormMLE(mleA,0)[1,]	  
    sdB <- 0.02
    gridx <- a <- seq(0,1,length.out=91)[-c(1,91)]
    gridy <- seq(-0.2,1.1,length.out=51)
    gridX <- expand.grid(gridx, gridy)
    #
    logden0 <- denBanana(c(mleA,0), coefA=coefA, sdB=sdB)
    luden <- apply( gridX, 1, denBanana, coefA=coefA, sdB=sdB)
    mLuden <- matrix(luden,nrow=length(gridx))
    imax <- which( matrix(luden,nrow=length(gridx))==max(luden), arr.ind=TRUE)[1,]
    
    
    blocks <- list(
            m1 = blockSpec(,,   # empty defaults to all parameters 
                    new("MetropolisBlockUpdater",
                            fLogDensity=denBanana
                            ,argsFLogDensity=list(coefA=coefA, sdB=sdB)
                            ,logDensityComponentNames = names(logden0)
                    )
            )
    )  
    theta0 <- c(a=0.8,b=0.8) 
    covarTheta0 <- diag(c(a=0.05,b=0.1)) 		
    sampler0 <- newPopulationSampler( blocks, theta0, covarTheta0, nPopulation=2L, nChainInPopulation=9L, nChainGroupsInPopulation=3L )
    sampler <- setupAndSample(sampler0, nSample=40L)
    
    sampler <- setupAndSample(sampler, nSample=80L)
    
    #plot( asMcmc(getSampleLogs(sampler)), smooth=FALSE)
    sampleLogsTail <- subsetTail(getSampleLogs(sampler),0.6)
    #plot( asMcmc(sampleLogsTail), smooth=FALSE )
    aC <- getAcceptanceTracker(sampler)
    computePopulationAcceptanceRates(aC)
    #
    stacked <- stackChainsInPopulation(sampleLogsTail)    
    ss <- ss0 <- t(getParametersForAllChains(sampleLogsTail))    
    
    image( gridx, gridy,  matrix(exp(luden),nrow=length(gridx)), col = rev(heat.colors(100)), xlab="a", ylab="b" )
    points( gridx[ imax[1] ], gridy[ imax[2] ], pch="X"  )
    points( ss[,1], ss[,2] )
    
    # here with no groupings
    sampler0UnGrouped <- newPopulationSampler( blocks, theta0, covarTheta0, nPopulation=2L, nChainInPopulation=3L )
    samplerUnGrouped <- setupAndSample(sampler0UnGrouped, nSample=40L*3L)
    #plot( asMcmc(getSampleLogs(samplerUnGrouped)), smooth=FALSE)
    
    samplerUnGrouped <- setupAndSample(samplerUnGrouped, nSample=80L*3L)
    
    sampleLogsTailUnGrouped <- subsetTail(getSampleLogs(samplerUnGrouped),0.6)
    #plot( asMcmc(sampleLogsTailUnGrouped), smooth=FALSE )
    aCUngrouped <- getAcceptanceTracker(samplerUnGrouped)
    computePopulationAcceptanceRates(aCUngrouped)
    
    .tmp.inspectAcceptanceRate <- function(){
        matplot(getAcceptanceLog(aC)[1L,,], type="l")
        matplot(getAcceptanceLog(aCUngrouped)[1L,,], type="l")
        hist(getAcceptanceLog(aC)[1L,,])
        hist(getAcceptanceLog(aCUngrouped)[1L,,])
    }
    
    #
    ssUngrouped <- ss0 <- t(getParametersForAllChains(sampleLogsTailUnGrouped))    
    
    image( gridx, gridy,  matrix(exp(luden),nrow=length(gridx)), col = rev(heat.colors(100)), xlab="a", ylab="b" )
    points( gridx[ imax[1] ], gridy[ imax[2] ], pch="X"  )
    points( ssUngrouped[,1], ssUngrouped[,2] )
    
    plot( density(ssUngrouped[,1]), xlim=c(0,1) )
    pX <- dlogitnorm(gridx, coefA[1], coefA[2])
    lines( pX ~ gridx, col="orange", type="l")
    lines( density(ss[,1]), col="blue")
    
    
    
}


