#require(testthat)
context("BatchSampler")

x0 <- c(a=3,b=5,c=7)
blocks = list(
        m1=blockSpecMock(c("a"),,c("obs1"))
        ,m2Parms =blockSpecMock(c("b","c"),,c("obs2","parms"))
)
sampler <- new("BatchSampler"
        , chainSampler=newChainSampler(blocks, names(x0)))
chainSampler <- getChainSampler(sampler)


sampleLogs(sampler) <- createSampleLogsWithInitialPopulation(sampler
        ,thetaPrior=x0
        ,covarTheta=diag( x0*0.1)
)


#chainSampler <- getChainSampler(sampler)
#sDim <- getSampleDimensions(sampler)
#nInterval <- 3L
#samplerMock <- initializeSampleDimensions(new("BatchSamplerMock", chainSampler=chainSampler )
#    , nPopulation = getNPopulation(sDim)
#    , nChainInPopulation=getNChainInPopulation(sDim)
#    , thin=getThin(sDim)
#    , nIntervalPerRange=nInterval
#)



test_that("newBatchSampler",{
            sampler <- newBatchSampler(blocks, x0, diag(0.1*x0), thin=2L) 
            nSampleBatchPopulations <- c(2L,3L)
            tmpDir <- tempdir()
            on.exit( unlink(tmpDir), TRUE)
            updatedSampler <- sampleBatches(sampler
                    , nBatchMax=2L
                    , nSampleBatch=8L
                    , TStart=10
                    , TEnd=1  
                    , restartFilename = file.path(tmpDir,"test_BatchSampler")
            )
            resDir <- dir(tmpDir)
            expect_true( all(c("test_BatchSampler_2.RData","test_BatchSampler_end.RData") %in% resDir ))
            getBatchLog(updatedSampler)
            sampleLogs <- getSampleLogs(updatedSampler)
            parms <- getParametersForPopulation(sampleLogs,1L)
            expect_that( dim(parms), equals(c(length(x0),(1L+2L*8L),getNChainInPopulation(sampleLogs))))
            # here the mock updater is used and initial x0 is not modified
            expect_true( all(parms[,,1] ==x0) )
            getEndTemperaturesPopulations(sampleLogs)
            #plot( getStepTemperaturesPopulations(sampleLogs)[[1L]][1L,] )
        })



