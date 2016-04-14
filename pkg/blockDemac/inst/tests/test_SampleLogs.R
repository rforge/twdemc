#require(testthat)
context("SampleLogs")

thin <- 4L
nSamplePopulations <- c(24L,32L)
nChainInPopulation <- 4L

sDim <- .fixtureSampleDimensions(nSample=nSamplePopulations, nChainInPopulation = nChainInPopulation)
fs <- .fixtureSampleLogs(sDim)

# only one component in logDensity
fs1D <- .fixtureSampleLogs(.fixtureSampleDimensions(blockSpecifications = list(
            add7 = blockSpecMock("a",,)
            ,met1 = blockSpecMock(c("b","c"),,"obs1")
        )), logDensityComponents0 = c(obs1=-10))

# no logDensity components
fs0D <- .fixtureSampleLogs(.fixtureSampleDimensions(blockSpecifications = list(
            add7 = blockSpecMock(c("a","b","c"),,)
        )), logDensityComponents0 = c())

fs0Sample <- .fixtureSampleLogs(.fixtureSampleDimensions(nSample=c(0L,0L)))
#getParametersAndInitialForPopulation(fs0Sample$sampleLogs,1L)

nParameter <- length(fs$x0)

sampleLogs <- new("SampleLogs", sampleDimensions=sDim)

m0 <- computeNRepresentativeRows( nParameter, nChainInPopulation )
initialPopulation <- array(1:nParameter, dim=c(nParameter, m0, getNChain(sDim)))
sampleLogsI <- new("SampleLogs", sampleDimensions=sDim, initialPopulation=initialPopulation)



test_that("getSampleDimensions",{
        #expect_that( getThin(sampleLogs), equals(thin))
        nBlock <- getNBlock(sDim)
        expect_that( getNBlock(sampleLogs), equals(nBlock))
        expect_that( getBlockNames(sampleLogs), equals(getBlockNames(sDim)))
        expect_that( getNChain(sampleLogs), equals(getNChain(sDim)))
        expect_that( getNPopulation(sampleLogs), equals(2))
    })

test_that("nSample->",{
        sl <- fs$sampleLogs
        nSampleBefore <- getNSamplePopulations(sl)
        nSamplePopulations(sl) <- 30L
        expect_that( getNSamplePopulations(sl), equals( rep(30L,getNPopulation(sl))))
        parms1 <- getParametersForPopulation(sl, 1L)
        expect_that( dim(parms1), equals(c(nParameter, 30L,getNChainInPopulation(sl) )))
        expect_true( !any(is.na(parms1[,1:nSampleBefore[1L],])) )
        expect_true( all(is.na(parms1[,-(1:nSampleBefore[1L]),])) )
        #
        sl <- fs0Sample$sampleLogs
        nSampleBefore <- getNSamplePopulations(sl)
        nSamplePopulations(sl) <- 30L
        expect_that( getNSamplePopulations(sl), equals( rep(30L,getNPopulation(sl))))
        parms1 <- getParametersForPopulation(sl, 1L)
        expect_that( dim(parms1), equals(c(nParameter, 30L,getNChainInPopulation(sl) )))
        expect_true( all(is.na(parms1[,-(1:nSampleBefore[1L]),])) )
        
    })

test_that("initialPopulation",{
        parms1 <- getParametersForPopulation(sampleLogs, 2L)
        expect_that( dim(parms1), equals(c(nParameter, nSamplePopulations[2], nChainInPopulation)))
        expect_true( all(is.na(parms1)) )
        parms1 <- getParametersAndInitialForPopulation(sampleLogs, 2L)
        expect_that( dim(parms1), equals(c(nParameter, nSamplePopulations[2], nChainInPopulation)))
        expect_true( all(is.na(parms1)) )
        #
        parms1 <- getParametersForPopulation(sampleLogsI, 2L)
        expect_that( dim(parms1), equals(c(nParameter, nSamplePopulations[2], nChainInPopulation)))
        expect_true( all(is.na(parms1)) )
        #sampleLogs <- initializeSampleLogs(new("SampleLogs", sampleDimensions=sDim, initialPopulation=numeric(0)))
        sampleLogs <- new("SampleLogs", sampleDimensions=sDim)
        parms1I <- getParametersAndInitialForPopulation(sampleLogsI, 2L)
        expect_that( dim(parms1I), equals(c(nParameter, m0+nSamplePopulations[2], nChainInPopulation)))
        x0Init <- fs$x0; x0Init[] <- as.numeric(seq_along(x0Init))
        expect_that( parms1I[,2,3], equals(x0Init) )            
        expect_true( all(!is.na(parms1I[,m0,])) )           
    })

test_that("setting and getting stepTemperatures",{
        stepTemperatures <- getStepTemperaturesPopulations(sampleLogs)
        expect_that( length(stepTemperatures), equals( getNPopulation(sampleLogs)))
        stepTemperaturesPop2 <- stepTemperatures[[2L]] 
        expect_that( dim(stepTemperaturesPop2), equals(c(getNLogDensityComponent(sampleLogs),getNSamplePopulations(sampleLogs)[2L])) )
        expect_true( all(stepTemperaturesPop2 == 1) )
        #
        stepTemperatures0 <- getStepTemperaturesPopulations(sampleLogs)
        newTemp <- lapply( seq_along(stepTemperatures0), function(iPop){
                stepTemperatures0[[iPop]]*(2*iPop) })
        stepTemperaturesPopulations(sampleLogs, rep(0L, getNPopulation(sampleLogs))) <- newTemp
        stepTemperatures <- getStepTemperaturesPopulations(sampleLogs)
        expect_that( length(stepTemperatures), equals( getNPopulation(sampleLogs)))
        #iPop <- 2L
        for( iPop in 1:getNPopulation(sampleLogs) ){
            stepTemperaturesPop <- stepTemperatures[[iPop]] 
            expect_that( dim(stepTemperaturesPop), equals(c(getNLogDensityComponent(sampleLogs),getNSamplePopulations(sampleLogs)[iPop])) )
            expect_true( all(stepTemperaturesPop == 2*iPop) )
        }
    })



test_that("record chain sampling result ",{
        nSample <- 4L
        blockNames <- getBlockNames(sDim)
        sampleLogPopulation <- initializeLog( new("SampleLogPopulation")
            ,nSample=32, parameterNames=names(fs$x0), logDensityComponentNames=names(fs$logDensityComponents0)
            ,blockNames=blockNames, nChain=nChainInPopulation )
        sampleLog1 <- new("SampleLogChainImpl", nSample=nSample
            , parameters=fs$x0Chains[,1:nSample]
            , logDensityComponents=fs$logDensityComponents0Chains[,1:nSample]
            #, countAcceptedInInterval=matrix(1:getNBlock(sDim), nrow=getNBlock(sDim), ncol=nSample
            #    , dimnames=list(blockNames,NULL))
        )
        sampleLog2 <- new("SampleLogChainImpl", nSample=nSample
            , parameters=2*fs$x0Chains[,1:nSample]
            , logDensityComponents=2*fs$logDensityComponents0Chains[,1:nSample]
            #, countAcceptedInInterval=2*matrix(1:getNBlock(sDim), nrow=getNBlock(sDim), ncol=nSample
            #    , dimnames=list(blockNames,NULL))
        )
        sampleLogList <- c(sampleLog1, sampleLog2, sampleLog2, sampleLog1)
        updatedSampleLogPopulation <- recordSampleLogs(sampleLogPopulation, sampleLogList, 4L)
        expect_that( getParameters(updatedSampleLogPopulation)[,4+(1:nSample),1], equals(getParameters(sampleLog1)))
        expect_that( getParameters(updatedSampleLogPopulation)[,4+(1:nSample),2], equals(getParameters(sampleLog2)))
        expect_that( getLogDensityComponents(updatedSampleLogPopulation)[,4+(1:nSample),1], equals(getLogDensityComponents(sampleLog1)))
        expect_that( getLogDensityComponents(updatedSampleLogPopulation)[,4+(1:nSample),2], equals(getLogDensityComponents(sampleLog2)))
        #expect_that( getCountAcceptedInInterval(updatedSampleLogPopulation)[,4+(1:nSample),1], equals(getCountAcceptedInInterval(sampleLog1)))
        #expect_that( getCountAcceptedInInterval(updatedSampleLogPopulation)[,4+(1:nSample),2], equals(getCountAcceptedInInterval(sampleLog2)))
    })


test_that("record sampling result",{
        iPopulationsStep <- 1:2
        sampleLogChains <- c(getSampleLogChainsForPopulations(fs$sampleLogs,1L),getSampleLogChainsForPopulations(fs$sampleLogs,2L))
        sampleLogChain <- sampleLogChains[[1]]        
        sampleLogChains <- lapply(sampleLogChains, function(sampleLogChain){
                nSample(sampleLogChain) <- 4L   # only record 4 steps instead of all 24
                sampleLogChain
            })
        nSample <- getNSample(sampleLogChains[[1]])
        updatedSampleLogs <- recordSampleLogChains(sampleLogsI, sampleLogChains, iPopulations=iPopulationsStep )
        # in creation of fs, x0 is multiplied by chainInPop and step
        expect_that( getParametersForPopulation(updatedSampleLogs, 1L)[,nSample,1], equals(fs$x0*nSample))
        expect_that( getParametersForPopulation(updatedSampleLogs, 1L)[,nSample,2], equals(2*fs$x0*nSample))
        expect_that( getParametersForPopulation(updatedSampleLogs, 2L)[,nSample,1], equals(fs$x0*nSample))
        expect_that( getLogDensityComponentsForPopulation(updatedSampleLogs, 2L)[,nSample,1], equals(fs$logDensityComponents0*nSample))
        #getCountAcceptedInIntervalForPopulation(sampleLogs, 2L)
        #getCountAcceptedInIntervalForPopulation(updatedSampleLogs, 2L)
        #expect_that( getCountAcceptedInIntervalForPopulation(updatedSampleLogs, 2L)[,nSample,1]
        #    , equals(fs$countAcceptedInInterval0))
        #
        iSample0 <- 8L  # record at different position
        updatedSampleLogs <- recordSampleLogChains(sampleLogsI, sampleLogChains, iPopulations=iPopulationsStep, iSample0=c(iSample0,iSample0) )
        expect_that( getParametersForPopulation(updatedSampleLogs, 1L)[,iSample0+nSample,1], equals(fs$x0*nSample))
        expect_that( getParametersForPopulation(updatedSampleLogs, 2L)[,iSample0+nSample,1], equals(fs$x0*nSample))
    })

test_that("stackChainsInPopulation SampleLogPopulation",{
        fs <- .fixtureSampleLogs()
        sDim <- getSampleDimensions(fs$sampleLogs)
        popLog <- fs$sampleLogs@populationLogs[[1]]
        expN <- getNSample(popLog)*getNChain(popLog)
        #mergeMethod <- "slice"
        lapply( c("stack","random","slice"), function(mergeMethod){
                res <- stackChains(popLog, mergeMethod)
                expect_that( getNSample(res), equals(expN) )
                expect_that( ncol(getParameters(res)), equals(expN) )
                expect_that( ncol(getLogDensityComponents(res)), equals(expN) )
                #expect_that( ncol(getCountAcceptedInInterval(res)), equals(expN) )
                expect_that( ncol(getStepTemperatures(res)), equals(expN) )
            })
    })

test_that("stackChainsInPopulation SampleLogs",{
        fs <- .fixtureSampleLogs()
        sDim <- getSampleDimensions(fs$sampleLogs)
        s2 <- stackChainsInPopulation(fs$sampleLogs)
        #s2 <- stackChainsInPopulation(fs$sampleLogs, "slice")
        expN <- getNSamplePopulations(sDim)*getNChainInPopulation(sDim)
        sDim2 <- getSampleDimensions(s2) 
        expect_that( getNSamplePopulations(sDim2), equals(expN))
        expect_that( getNChainInPopulation(sDim2), equals(1L))
        expect_that( getNChain(sDim2), equals(getNPopulation(sDim2)))
        parameters <- getParametersForPopulation(s2,1L)
        expect_that( dim(parameters), equals( c(getNParameter(sDim), expN[1], 1L)))
        .tmp.f <- function(){
            plot(asMcmc(s2), smooth=FALSE)
        }
    })

test_that("subset SampleLogPopulation",{
        fs <- .fixtureSampleLogs()
        sDim <- getSampleDimensions(fs$sampleLogs)
        popLog <- fs$sampleLogs@populationLogs[[1]]
        iSample <- seq(1, min(getNSamplePopulations(sDim)), by=5)
        expN <- length(iSample)
        popLog2 <- subSample(popLog, iSample)
        expect_that( getNSample(popLog2), equals(expN) )
        expect_that( ncol(getParameters(popLog2)), equals(expN) )
        expect_that( ncol(getLogDensityComponents(popLog2)), equals(expN) )
        #expect_that( ncol(getCountAcceptedInInterval(popLog2)), equals(expN) )
        expect_that( ncol(getStepTemperatures(popLog2)), equals(expN) )
    })

test_that("subset SampleLogs",{
        fs <- .fixtureSampleLogs()
        sDim <- getSampleDimensions(fs$sampleLogs)
        iSample <- seq(1, min(getNSamplePopulations(sDim)), by=5)
        expN <- length(iSample)
        s2 <- subSample(fs$sampleLogs, iSample)
        sDim2 <- getSampleDimensions(s2) 
        expect_that( getNSamplePopulations(sDim2), equals(rep(expN, getNPopulation(sDim))))
        parameters <- getParametersForPopulation(s2,1L)
        expect_that( dim(parameters), equals( c(getNParameter(sDim), expN, getNChainInPopulation(sDim))))
        .tmp.f <- function(){
            plot(asMcmc(s2), smooth=FALSE)
        }
    })

test_that("squeeze SampleLogs",{
            fs <- .fixtureSampleLogs()
            sDim <- getSampleDimensions(fs$sampleLogs)
            nOrig <- getNSamplePopulations(sDim)
            fraction <- 0.5
            s2 <- squeeze(fs$sampleLogs, fraction=fraction )
            expN <- rep( round(min(nOrig) * fraction), getNPopulation(fs$sampleLogs))
            sDim2 <- getSampleDimensions(s2) 
            expect_that( getNSamplePopulations(sDim2), equals(expN))
            #expect_that( getThin(sDim2), equals(getThin(sDim)*2L))
            parameters <- getParametersForPopulation(s2,1L)
            expect_that( dim(parameters), equals( c(getNParameter(sDim), expN[1], getNChainInPopulation(sDim))))
            .tmp.f <- function(){
                plot(asMcmc(s2), smooth=FALSE)
            }
        })

test_that("squeeze SampleLogs, error on few samples",{
            fs <- .fixtureSampleLogs(nInitialSample=20L)
            sDim <- getSampleDimensions(fs$sampleLogs)
            nOrig <- getNSamplePopulations(sDim)
            expect_error(
                s2 <- squeeze(fs$sampleLogs, nSample=30 )
            )
        })


test_that("subsetTail SampleLogs",{
        fs <- .fixtureSampleLogs()
        sDim <- getSampleDimensions(fs$sampleLogs)
        s2 <- subsetTail(fs$sampleLogs, 0.8 )
        expN <- ceiling(getNSamplePopulations(sDim) *0.8)
        sDim2 <- getSampleDimensions(s2) 
        expect_that( getNSamplePopulations(sDim2), equals(expN))
        parameters <- getParametersForPopulation(s2,1L)
        expect_that( dim(parameters), equals( c(getNParameter(sDim), expN[1], getNChainInPopulation(sDim))))
        .tmp.f <- function(){
            plot(asMcmc(s2), smooth=FALSE)
        }
    })

test_that("subsetTail SampleLogs nSample argument",{
            fs <- .fixtureSampleLogs()
            sDim <- getSampleDimensions(fs$sampleLogs)
            nSamplePops <- getNSamplePopulations(sDim) 
            expect_warning(s2 <- subsetTail(fs$sampleLogs, nSample=28L ))
            expN <- pmin(nSamplePops, 28L)
            sDim2 <- getSampleDimensions(s2) 
            expect_that( getNSamplePopulations(sDim2), equals(expN))
            parameters <- getParametersForPopulation(s2,1L)
            expect_that( dim(parameters), equals( c(getNParameter(sDim), expN[1], getNChainInPopulation(sDim))))
            .tmp.f <- function(){
                plot(asMcmc(s2), smooth=FALSE)
            }
        })

test_that("asMcmc SampleLogs",{
            fs <- .fixtureSampleLogs()
            sDim <- getSampleDimensions(fs$sampleLogs)
            mc <- asMcmc(fs$sampleLogs)
            expect_true( inherits(mc,"mcmc.list"))
            expect_that( length(mc), equals(getNChain(sDim))) 
            expect_that( nrow(mc[[1]]), equals( min(getNSamplePopulations(sDim))))
            .tmp.f <- function(){
                plot(mc, smooth=FALSE)
            }
        })



test_that("append SampleLogPopulation",{
        fs <- .fixtureSampleLogs()
        sDim <- getSampleDimensions(fs$sampleLogs)
        sl1 <- fs$sampleLogs@populationLogs[[1]]
        sl2 <- fs$sampleLogs@populationLogs[[2]]
        slc <- appendSampleLog(sl1, sl2)
        expect_that( getNSample(slc), equals(getNSample(sl1)+getNSample(sl2)))
        expect_that( getParameters(slc)[,1:getNSample(sl1),], equals(getParameters(sl1)))
        expect_that( getParameters(slc)[,getNSample(sl1)+1:getNSample(sl2),], equals(getParameters(sl2)))
        .tmp.f <- function(){
            plot(asMcmc(s2), smooth=FALSE)
        }
    })

test_that("setting nSample SampleLogPopulation",{
        fs <- .fixtureSampleLogs()
        sDim <- getSampleDimensions(fs$sampleLogs)
        sl1 <- fs$sampleLogs@populationLogs[[1]]
        nSampleBefore <- getNSample(sl1)
        nSample <- nSampleBefore + 5L
        nSample(sl1) <- nSample
        expect_that( getNSample(sl1), equals(nSample))
        nSample <- nSampleBefore - 2L
        nSample(sl1) <- nSample
        expect_that( getNSample(sl1), equals(nSample))
    })


test_that("setting nSample SampleLogs",{
        fs <- .fixtureSampleLogs()
        sDim <- getSampleDimensions(fs$sampleLogs)
        sl <- fs$sampleLogs
        #getNSamplePopulations(fs$sampleLogs)
        nSample <- 32L
        nSamplePopulations(sl) <- nSample  
        expect_that(getNSamplePopulations(sl), equals(rep(nSample,getNPopulation(sl))))
        expect_that(ncol(getLogDensityComponentsForPopulation(sl,1L)), equals(nSample))
    })


test_that("applyPopulations SampleLogs",{
        fs <- .fixtureSampleLogs()
        parameters <- applyPopulations(fs$sampleLogs, getParametersForPopulation)
        expect_true( is.list(parameters) )
        expect_that( length(parameters), equals(getNPopulation(fs$sampleLogs)) )
    })

test_that("computeTemperatedLogDensity SampleLogs",{
        fs <- .fixtureSampleLogs()
        sDim <- getSampleDimensions(fs$sampleLogs)
        logDenT <- computeTemperatedLogDensityForPopulation(fs$sampleLogs, 2L)
        expect_that( dim(logDenT), equals(dim(getLogDensityComponentsForPopulation(fs$sampleLogs,2L))))        
    })

test_that("computeBlockLogDensities SampleLogsPopulation",{
            fs <- .fixtureSampleLogs()
            sampleLogs <- fs$sampleLogs
            sDim <- getSampleDimensions(sampleLogs)
            logPop <- getSampleLogOfPopulation(sampleLogs, 2L)
            logDenBlocks <- computeBlockLogDensities( logPop, sampleLogs )
            expect_that( dim(logDenBlocks), equals(c(2L,getNSamplePopulations(sampleLogs)[2L],getNChainInPopulation(sampleLogs))))
            expect_that( rownames(logDenBlocks), equals(getBlockNames(sampleLogs)) )
            logDenT <- computeTemperatedLogDensityForPopulation(sampleLogs, 2L)
            exp <- apply(logDenT[2:3,,], 2:3, sum); dimnames(exp) <- NULL
            expect_that( logDenBlocks[2L,,], equals(exp) )
            #
            # only second block uses logDensityComponents
            sampleLogs <- fs1D$sampleLogs
            logPop <- getSampleLogOfPopulation(sampleLogs, 2L)
            logDenBlocks <- computeBlockLogDensities( logPop, sampleLogs )
            expect_that( dim(logDenBlocks), equals(c(1L,getNSamplePopulations(sampleLogs)[2L],getNChainInPopulation(sampleLogs))))
            expect_that( rownames(logDenBlocks), equals(getBlockNames(sampleLogs)[2L]) )
            #
            # none of the blocks uses logDensityComponents
            sampleLogs <- fs0D$sampleLogs
            logPop <- getSampleLogOfPopulation(sampleLogs, 2L)
            expect_warning(
                logDenBlocks <- computeBlockLogDensities( logPop, sampleLogs )
            )
            expect_that( dim(logDenBlocks), equals(c(0L,getNSamplePopulations(sampleLogs)[2L],getNChainInPopulation(sampleLogs))))
            expect_that( length(rownames(logDenBlocks)), equals(0L) )
        })


test_that("computeSampleRanksForPopulation SampleLogs",{
        fs <- .fixtureSampleLogs()
        sampleRanks <- computeSampleRanksForPopulation(fs$sampleLogs,2L)
        expect_that( dim(sampleRanks), equals(c(getNSamplePopulations(fs$sampleLogs)[2L],getNChainInPopulation(fs$sampleLogs))))
        expect_that( sampleRanks[,1], equals(1:getNSamplePopulations(fs$sampleLogs)[2L]))
        #
        sampleLogs <- fs1D$sampleLogs
        sampleRanks <- computeSampleRanksForPopulation(sampleLogs,2L)
        expect_that( dim(sampleRanks), equals(c(getNSamplePopulations(sampleLogs)[2L],getNChainInPopulation(sampleLogs))))
        expect_that( sampleRanks[,1], equals(1:getNSamplePopulations(sampleLogs)[2L]))
        #
        sampleLogs <- fs0D$sampleLogs
        expect_warning( sampleRanks <- computeSampleRanksForPopulation(sampleLogs,2L) )
        expect_that( dim(sampleRanks), equals(c(getNSamplePopulations(sampleLogs)[2L],getNChainInPopulation(sampleLogs))))
        # here default order should be returned
        expect_that( sampleRanks[,1], equals(1:getNSamplePopulations(sampleLogs)[2L]))
    })

test_that("getBestSamples SampleLogs",{
            popLog <- fs$sampleLogs@populationLogs[[1]]
            popLogBest <- getBestSamples(popLog, fs$sampleDimensions, nBest=3L)
            logDenComp <- getLogDensityComponents(popLogBest)[,,1]
            #getLogDensityComponents(popLog)[,1:2,1]
            # here the best ranks are two first stpes of the first chain, and the first of the second chain
            expect_that( logDenComp, equals(cbind(
                  getLogDensityComponents(popLog)[,1:2,1]
                    ,getLogDensityComponents(popLog)[,1,2]
            )))
            #
            popLog <- fs1D$sampleLogs@populationLogs[[1]]
            popLogBest <- getBestSamples(popLog, fs1D$sampleDimensions, nBest=3L)
            logDenComp <- adrop(getLogDensityComponents(popLogBest)[,,1 ,drop=FALSE],3)
            expect_that( logDenComp, equals(cbind(
                                    adrop( getLogDensityComponents(popLog)[,1:2,1, drop=FALSE],3)
                                    ,adrop(getLogDensityComponents(popLog)[,1,2, drop=FALSE],3)
                            )))
            #
            popLog <- fs0D$sampleLogs@populationLogs[[1]]
            expect_warning( popLogBest <- getBestSamples(popLog, fs0D$sampleDimensions, nBest=3L) )
            parms <- getParameters(popLogBest)
            expect_that( parms, equals(getParameters(popLog)[,1:3,1 ,drop=FALSE]))
        })






