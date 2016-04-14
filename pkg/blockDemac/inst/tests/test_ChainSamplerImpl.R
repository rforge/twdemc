#require(testthat)
context("ChainSamplerImpl")

# Fixture is reused by different test cases
fx <- .fixtureChainSampler()
blockUpdaters <- newBlockUpdaters(fx$blocks, list(), names(fx$x0))
chainSampler <- new("ChainSamplerImpl", blockUpdaters=blockUpdaters, thin=fx$thin)
chainSampler <- setRangeSpecs(chainSampler,
        chainState=fx$chainState, nInterval=fx$nInterval, stepInfoRange=fx$stepInfoRange)


test_that("initialization",{rt
            chainSampler0 <- new("ChainSamplerImpl")
            chainSampler <- initializeBlockUpdaters(chainSampler0, blockUpdaters, fx$thin)
                chainSampler <- initializeBlockUpdaters(chainSampler0, blockUpdaters) # missing thin
            expect_error(            
                    chainSampler <- initializeBlockUpdaters(chainSampler0, blockUpdaters, thin=c(1L,1L) ) # missing thin
            )
        })


test_that(".computeUpdatedChainStateThinningIntervalAndCountAccepted",{
            #getChainState(chainSampler)
            updatedChainSampler <- .computeUpdatedChainStateThinningIntervalAndCountAccepted(chainSampler)
            x1 <- getChainStateParameters(getChainState(updatedChainSampler))
            expect_that( x1, equals(fx$x0 +c(a=7*fx$thin, isM=0, m=7*fx$thin, isC=0, c=9*fx$thin)))
            cntAccepted <- updatedChainSampler@countAcceptedInCurrentInterval
            expect_that( cntAccepted, equals(rep(fx$thin,3)))
            logDensityComponents <- getLogDensityComponents(getChainState(updatedChainSampler))
            # m1 invalidated by following block
            expect_that( logDensityComponents, equals(c(m1=NA,m2=-10,parms=-2)))
        })

test_that(".sampleInterval",{
            updatedChainSampler <- .sampleInterval(chainSampler)
            #getChainState(updatedChainSampler)
            sampleLog <- getSampleLog(updatedChainSampler)
            parms <- getParameters(sampleLog)
            expect_that( parms[,1], equals(fx$x0 +c(a=7*fx$thin, isM=0, m=7*fx$thin, isC=0, c=9*fx$thin)))
            logDensityComponents <- getLogDensityComponents(sampleLog)
            expect_that( logDensityComponents[,1], equals(c(m1=-10,m2=-10,parms=-2)))
            proportionAcceptedInInterval <- getProportionAcceptedInInterval(sampleLog)
            expect_that( proportionAcceptedInInterval[,1], is_equivalent_to(rep(fx$thin/fx$thin,3)))
        })

test_that("sampleRange",{
            updatedChainSampler <- sampleRange(chainSampler)
            #getChainState(updatedChainSampler)
            sampleLog <- getSampleLog(updatedChainSampler)
            parms <- getParameters(sampleLog)
            expect_that( parms[,3], equals(fx$x0 +3*c(a=7*fx$thin, isM=0, m=7*fx$thin, isC=0, c=9*fx$thin)))
            proportionAcceptedInInterval <- getProportionAcceptedInInterval(sampleLog)
            expect_that( proportionAcceptedInInterval[,3], is_equivalent_to(rep(fx$thin/fx$thin,3)))
            logDensityComponents <- getLogDensityComponents(sampleLog)
            expect_that( logDensityComponents[,3], equals(c(m1=-10,m2=-10,parms=-2)))
        })

test_that("setting subSpace",{
            blockUpdater1 <- getBlockUpdaterForIndex(chainSampler@blockUpdaters, 1L)
            resSubSpace <- getSubSpace(blockUpdater1)
            expect_that( getLowerParBounds(resSubSpace), equals(numeric(0)))
            # indexedBounds initialized when creating BlockUpdater
            parameterNames <- getParameterNames(getBlockDimensions(chainSampler))
            expect_that( rownames(getIndexedParBounds(resSubSpace)), equals(parameterNames))
            #
            subSpaceAll <- new("SubSpace")
            parameterNames <- getParameterNames(getBlockDimensions(chainSampler))
            subSpaceAll <- initializeIndexedBoundsForParameterNames(subSpaceAll, parameterNames)
            expect_that( rownames(getIndexedParBounds(subSpaceAll)), equals(parameterNames))
            subSpaceUpper <- getSplitSubSpaces(subSpaceAll, c(a=10) )$upper
            expect_that( rownames(getIndexedParBounds(subSpaceUpper)), equals(parameterNames))
            subSpace(chainSampler) <- subSpaceUpper            
            blockUpdater1 <- getBlockUpdaterForIndex(chainSampler@blockUpdaters, 1L)
            resSubSpace <- getSubSpace(blockUpdater1)
            expect_that( getLowerParBounds(resSubSpace), equals(getLowerParBounds(subSpaceUpper)))
            expect_that( getIndexedParBounds(subSpaceUpper)["a","lower"], equals(10))
        })
    
test_that("sampling linear model",{
        fxLinReg <- .fixtureLinReg1()
        lmDummy <- with(fxLinReg, lm( obs ~ xval, weights=1/sdObs^2))		# results without priors
        (.expTheta <- coef(lmDummy))
        (.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
        (.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=names(fxLinReg$theta0)) )
        #
        blockUpdaters <- newBlockUpdaters(fxLinReg$blocks, list(), names(fxLinReg$theta0))
        thin <- 2L
        chainSampler <- new("ChainSamplerImpl", blockUpdaters=blockUpdaters, thin=thin)
        #        
        nGen <- 150L
        steps <- t(rmvnorm( nGen, mean=c(a=0,b=0), sigma=.expCovTheta*1.2 ))
        stepInfoRange <- new("StepInfoRangeImpl", step=steps)
        # starting conditions not in confidence interval 
        chainState <- new("ChainState", parameters=c(a=80,b=-20), logDensityComponents=c(obs1=NA_real_))
        #chainState <- new("ChainState", parameters=fxLinReg$theta0, logDensityComponents=c(obs1=NA_real_))
        chainSampler <- setRangeSpecs(chainSampler,
            chainState=chainState, nInterval=nGen%/%thin, stepInfoRange=stepInfoRange)
        #
        updatedChainSampler <- sampleRange(chainSampler)
        sampleLog <- getSampleLog(updatedChainSampler)
        parms <- t(getParameters(sampleLog))
        parmsBurnedIn <- parms[-(1:40),]
        .tmp.f <- function(){
            matplot(parms, type="l")
            matplot(parmsBurnedIn, type="l")
        }
        .mean <- colMeans(parmsBurnedIn)
        .sd <- apply(parmsBurnedIn, 2, sd)
        expect_isOfMagnitude( .expSdTheta, .sd )
        .pthetaTrue <- pnorm(.expTheta, mean=.mean, sd=.sd)
        expect_isInInterval( .pthetaTrue ) 
    })   



