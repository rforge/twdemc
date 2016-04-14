#require(testthat)
context("PopulationSampler")

x0 <- c(a=3,b=5,c=7)
blocks = list(
        m1=blockSpecMock(c("a"),,c("obs1"))
        ,m2Parms =blockSpecMock(c("b","c"),,c("obs2","parms"))
)
sampler <- new("PopulationSampler"
        , chainSampler=newChainSampler(blocks, names(x0)))
chainSampler <- getChainSampler(sampler)


sampleLogs(sampler) <- createSampleLogsWithInitialPopulation(sampler
        ,thetaPrior=x0
        ,covarTheta=diag( x0*0.1)
)


#chainSampler <- getChainSampler(sampler)
#sDim <- getSampleDimensions(sampler)
#nInterval <- 3L
#samplerMock <- initializeSampleDimensions(new("PopulationSamplerMock", chainSampler=chainSampler )
#    , nPopulation = getNPopulation(sDim)
#    , nChainInPopulation=getNChainInPopulation(sDim)
#    , thin=getThin(sDim)
#    , nIntervalPerRange=nInterval
#)



.tmp.f <- function(){
    populationSampler <- new("PopulationSamplerMock"
            ,populationSpecs = populationSpecs
            ,chainSampler = fs$chainSampler
            ,jumpProposer = fs$jumpProposer
            ,sampleLogs = fs$sampleLogs 
    )
}


test_that("createSampleLogsWithInitialPopulation",{
            nChainInPopulation <- 4L
            nPopulation <- 3L
            sampler <- new("PopulationSampler", chainSampler=chainSampler)
            sampleLogs <- createSampleLogsWithInitialPopulation(sampler
                    ,thetaPrior=x0
                    ,covarTheta=diag( x0*0.1)
                    ,nChainInPopulation = nChainInPopulation
                    ,nPopulation = nPopulation
            )
            expect_that( getNSamplePopulations(sampleLogs), equals(rep(1L,nPopulation)) )
            expect_that( getParameterNames(sampleLogs), equals(names(x0)) )
            parms1 <- getParametersForPopulation(sampleLogs,1L)
            expect_that(dim(parms1), equals(c(length(x0), 1L, nChainInPopulation)))
            expect_that(parms1[,1,1], equals(x0))
            parmsAndInitial1 <- getParametersAndInitialForPopulation(sampleLogs,1L)
            m0 <- computeNRepresentativeRows( length(x0),nChainInPopulation)
            expect_true( ncol(parmsAndInitial1) >= m0 )
            #
            logDen1 <- getLogDensityComponentsForPopulation(sampleLogs,1L)
            expect_that(dim(logDen1), equals(c(getNLogDensityComponent(sampleLogs), 1L, nChainInPopulation)))
            for( iPop in 1:getNPopulation(sampleLogs)){
                expect_true( all(is.finite(getLogDensityComponentsForPopulation(sampleLogs,iPop))) )
            }
        })

test_that("createSampleLogsWithInitialPopulation with nonfinite density components",{
            blocks = list(
                    m1=blockSpec(,,new("MetropolisBlockUpdater"
                                    , fLogDensity=function(parms, intermediates, logDen){
                                        c(obs1 = if( parms[1] < x0[1] ) -Inf else -10 ) 
                                    }
                                    , logDensityComponentNames = c("obs1")
                            ))
            )
            sDim <- .fixtureSampleDimensions(blockSpecifications=blocks)
            chainSampler <- newChainSampler(blocks, names(x0))
            samplerMock <- new("PopulationSamplerMock", chainSampler=chainSampler )
            x0 <- c(a=1,b=100,c=1e5)
            # may get to few finite sample by chance, try twice
            initialPopulation <- try( createSampleLogsWithInitialPopulation(samplerMock
                    ,thetaPrior=x0
                    ,covarTheta=diag( 0.2*x0 )
                    ,nChainInPopulation=getNChainInPopulation(sDim)
            ), silent=TRUE)
            if( inherits(initialPopulation, "try-error") )
                initialPopulation <-  createSampleLogsWithInitialPopulation(samplerMock
                                ,thetaPrior=x0
                                ,covarTheta=diag( 0.2*x0 )
                                ,nChainInPopulation=getNChainInPopulation(sDim)
                        )
            sampleLogs(samplerMock) <- initialPopulation 
            sampleLogs <- getSampleLogs(samplerMock)
            parms1 <- getParametersForPopulation(sampleLogs,1L)
            parms1i <- getParametersAndInitialForPopulation(sampleLogs,1L)
            for( iPop in 1:getNPopulation(sampleLogs)){
                expect_true( all(is.finite(getLogDensityComponentsForPopulation(sampleLogs,iPop))) )
            }
        })



#test_that("initializeSampleDimensions",{
#        #.samplePopulations(populationSampler)
#        sampler <- initializeSampleDimensions(new("PopulationSampler", chainSampler=chainSampler))
#        sampleLogs <- getSampleLogs(sampler)
#        expect_true( getNPopulation(sampleLogs) > 0 )
#        expect_true( all(getNSamplePopulations(sampleLogs) > 0) )
#        expect_true( getNChainInPopulation(sampleLogs) > 0 )
#        expect_that( getNChain(sampleLogs), equals(getNPopulation(sampleLogs)*getNChainInPopulation(sampleLogs)) )
#        expect_that( rownames(getParametersForPopulation(sampleLogs,1L)), equals(names(x0)) )
#        expect_that( length(getSubSpacePopulations(sampler)), equals(getNPopulation(sampleLogs)))
#        expect_that( length(sampler@stepTemperaturesPopulations), equals(getNPopulation(sampleLogs)))
#        #
#        sampler <- initializeSampleDimensions(sampler
#            ,parameterNames = names(x0)
#            ,blocks = blocks
#            , nSamplePopulations=rep(32L,2L)
#            , nChainInPopulation=3L
#            , thin=4L
#            , nIntervalPerRange=2L
#        )
#        sampleLogs <- getSampleLogs(sampler)
#        expect_that( getNSamplePopulations(sampleLogs), equals(rep(32L,2L)) )
#        expect_that( getNPopulation(sampleLogs), equals(2L) )
#        expect_that( getNChainInPopulation(sampleLogs), equals(3L) )
#        expect_that( getThin(sampleLogs), equals(4L) )
#        expect_that( getNChain(sampleLogs), equals(getNPopulation(sampleLogs)*getNChainInPopulation(sampleLogs)) )
#        #
#        expect_that( sampler@chainSampler@thin, equals(4L) )
#        expect_true( sampler@jumpProposer@deSettings["nRepresentativeRows"] > 0 )
#    })


test_that("set initialial population",{
            sampler <- sampler <- new("PopulationSampler", chainSampler=chainSampler)
            sampleLogs <- createSampleLogsWithInitialPopulation(sampler
                    ,thetaPrior=x0
                    ,covarTheta=diag( x0*0.1)
            )
            expect_true( !any(is.na(getStepTemperaturesPopulations(sampleLogs)[[2L]])) )
            sampleLogs(sampler) <- sampleLogs
        })


test_that(".initializeStepTemperatures",{
            nDenComp <- getNLogDensityComponent(getSampleDimensions(sampler))
#        sampleLogs <- getSampleLogs(sampler)
#        logDenNamed <- function(x){ structure(x,names=getLogDensityComponentNames(sampleLogs))}
#        expect_that( length(sampler@stepTemperaturesPopulations), equals(getNPopulation(sampleLogs)))
#        expect_that( getStepTemperaturesPopulations(sampleLogs)[[2L]][,1], equals(logDenNamed(rep(1,nDenComp))))
#        expect_that( getEndTemperaturesPopulations(sampleLogs)[,2], equals(logDenNamed(rep(1,nDenComp))))
            sampler <- sampler <- new("PopulationSampler", chainSampler=chainSampler)
            sampleLogs(sampler) <- createSampleLogsWithInitialPopulation(sampler
                    ,thetaPrior=x0
                    ,covarTheta=diag( x0*0.1)
            )
            sampler@nSampleBeforeSampling <- getNSamplePopulations(sampler@sampleLogs)
            sampler@nSampleBatchPopulations <- c(26L,32L)         
            nSamplePopulations(sampler@sampleLogs) <- getNSamplePopulations(sampler@sampleLogs) +  sampler@nSampleBatchPopulations
            updatedSampler <- .initializeStepTemperatures(sampler
                    ,TStart = 1+(1:nDenComp)
                    ,TEnd = rep(2, nDenComp)
            )
            sampleLogs <- getSampleLogs(updatedSampler)
            logDenNamed <- function(x){ structure(x,names=getLogDensityComponentNames(sampleLogs))}
            # use tolerance, because first step is already first decay step below specified TStart value
            expect_that( dim(getStepTemperaturesPopulations(sampleLogs)[[2L]]), equals(c(getNParameter(sampleLogs),getNSamplePopulations(sampleLogs)[2L])) )
            expect_true( !any(is.na(getStepTemperaturesPopulations(sampleLogs)[[2L]])) )
            expect_equal( getStepTemperaturesPopulations(sampleLogs)[[2L]][,sampler@nSampleBeforeSampling[2L]+1L], logDenNamed(1+(1:nDenComp)), tolerance = .2 )
            expect_that( getEndTemperaturesPopulations(sampleLogs)[,2], equals(logDenNamed(rep(2,nDenComp))))
        })



test_that("setupAndSample",{
            sampler <- sampler <- new("PopulationSampler", chainSampler=chainSampler)
            sampleLogs <- createSampleLogsWithInitialPopulation(sampler
                    ,thetaPrior=x0
                    ,covarTheta=diag( x0*0.1)
            )
            sampleLogs(sampler) <- sampleLogs
            nSampleBatchPopulations <- c(2L,3L)
            nDenComp <- getNLogDensityComponent(getSampleDimensions(sampler))
            TS <- 1+(1:nDenComp)
            updatedSampler <- setupAndSample(sampler, nSampleBatchPopulations, TStart = TS, thin=2L)
            #
            sampleLogs <- getSampleLogs(updatedSampler)
            parms <- getParametersForPopulation(sampleLogs,2L)
            expect_that( dim(parms), equals(c(getNParameter(sampleLogs),getNSamplePopulations(sampleLogs)[2L],getNChainInPopulation(sampleLogs))))
            expect_true( all(is.finite(parms)) )
            pAcceptPops <- getAcceptanceLog(getAcceptanceTracker(updatedSampler))
            expect_that( dim(pAcceptPops), equals(c(getNBlock(sampleLogs),max(nSampleBatchPopulations),getNChain(sampleLogs))))
            expect_true( all(is.finite(pAcceptPops[,min(nSampleBatchPopulations),])) )
            expect_true( all(is.na(pAcceptPops[,-(1:min(nSampleBatchPopulations)),1:2])) )
            expect_that( getEndTemperaturesPopulations(sampleLogs), equals(
                            structure(matrix(TS, nrow=length(TS), dimnames=list(getLogDensityComponentNames(sampleLogs),NULL), ncol=2L))
                    ))
            .tmp.f <- function(){
                matplot(1:ncol(parms), t(parms[,,1]), type="l")
                matplot(1:ncol(pAcceptPops), t(pAcceptPops[,,1]), type="l")
            }
        })


test_that("setupAndSample noBurnin",{
            sampler <- sampler <- new("PopulationSampler", chainSampler=chainSampler)
            sampleLogs <- createSampleLogsWithInitialPopulation(sampler
                    ,thetaPrior=x0
                    ,covarTheta=diag( x0*0.1)
            )
            sampleLogs(sampler) <- sampleLogs
            nSampleBatchPopulations <- c(1L,1L)
            isBurnin(sampler) <- FALSE
            updatedSampler <- setupAndSample(sampler, nSampleBatchPopulations)
            #
            nPop <- length(nSampleBatchPopulations)
            j2 <- updatedSampler@jumpProposer             
            #expect_that( .getNStepBackPops(j2, rep(0.25, nPop)), equals(rep(.Machine$integer.max,nPop)) )
            expect_that(isBurnin(j2), equals(FALSE) )
        })


test_that("setupAndSample: no Metropolis updaters (no logDen, no jumps)",{
            blocks = list(
                    m1=blockSpecMock(c("a"),,)
                    ,m2Parms =blockSpecMock(c("b","c"),,)
            )
            sDim <- .fixtureSampleDimensions(blockSpecifications=blocks)
            fs <- .fixtureSampleLogs(sampleDimensions=sDim, logDensityComponents0=numeric(0) )
            nInterval <- 3L
            chainSampler <- newChainSampler(blocks)
            #
            samplerMock <- new("PopulationSamplerMock", chainSampler=chainSampler, nIntervalPerRange=nInterval )
            sampleLogs(samplerMock) <- initialPopulation <- createSampleLogsWithInitialPopulation(samplerMock
                    ,thetaPrior=x0
                    ,covarTheta=diag( x0*0.1)
                    ,nChainInPopulation=getNChainInPopulation(sDim)
            )  
            sampleLogs <- getSampleLogs(samplerMock)
            for( iPop in 1:getNPopulation(sampleLogs)){
                expect_true( all(is.finite(getLogDensityComponentsForPopulation(sampleLogs,iPop) )) )
            }
            expect_warning(updatedSampler <- setupAndSample(samplerMock
                    , nSample=getNSamplePopulations(sDim)
            ))
            #
            sampleLogs <- getSampleLogs(updatedSampler)
            parms <- getParametersForPopulation(sampleLogs,2L)
            expect_that( dim(parms), equals(c(getNParameter(sampleLogs),getNSamplePopulations(sampleLogs)[2L],getNChainInPopulation(sampleLogs))))
            expect_true( all(is.finite(parms)) )
            pAcceptPops <- getAcceptanceLog(getAcceptanceTracker(updatedSampler))
            expect_that( dim(pAcceptPops), equals(c(getNBlock(sampleLogs),max(getNSamplePopulations(sDim)),getNChain(sampleLogs))))
            expect_true( all(is.finite(pAcceptPops[,min(getNSamplePopulations(sDim)),])) )
            expect_true( all(is.na(pAcceptPops[,-(1:min(getNSamplePopulations(sDim))),1:2])) )
        })

test_that("setupAndSample: 1D logDen, 1D parm",{
            blocks = list(
                    m1=blockSpecMock(c("a"),,c("obs1"))
            )
            sDim <- .fixtureSampleDimensions(blockSpecifications=blocks)
            fs <- .fixtureSampleLogs(sampleDimensions=sDim, logDensityComponents0=c(obs1=-10), x0=c(a=2) )
            nInterval <- 3L
            chainSampler <- newChainSampler(blocks)
            #
            samplerMock <- new("PopulationSamplerMock", chainSampler=chainSampler, nIntervalPerRange=nInterval )
            sampleLogs(samplerMock) <- initialPopulation <- createSampleLogsWithInitialPopulation(samplerMock
                    ,thetaPrior=fs$x0
                    ,covarTheta=diag( fs$x0*0.2, nrow=1)
                    ,nChainInPopulation=getNChainInPopulation(sDim)
            )  
            sampleLogs <- getSampleLogs(samplerMock)
            for( iPop in 1:getNPopulation(sampleLogs)){
                expect_true( all(is.finite(getLogDensityComponentsForPopulation(sampleLogs,iPop) )) )
            }
            updatedSampler <- suppressWarnings(
                    updatedSampler <- setupAndSample(samplerMock
                            , nSample=getNSamplePopulations(sDim))
                    )
            #
            sampleLogs <- getSampleLogs(updatedSampler)
            parms <- getParametersForPopulation(sampleLogs,2L)
            expect_that( dim(parms), equals(c(getNParameter(sampleLogs),getNSamplePopulations(sampleLogs)[2L],getNChainInPopulation(sampleLogs))))
            expect_true( all(is.finite(parms)) )
            pAcceptPops <- getAcceptanceLog(getAcceptanceTracker(updatedSampler))
            expect_that( dim(pAcceptPops), equals(c(getNBlock(sampleLogs),max(getNSamplePopulations(sDim)),getNChain(sampleLogs))))
            expect_true( all(is.finite(pAcceptPops[,min(getNSamplePopulations(sDim)),])) )
            expect_true( all(is.na(pAcceptPops[,-(1:min(getNSamplePopulations(sDim))),1:2])) )
        })

test_that("setupAndSample: extending existing sampleLogs",{
            fs <- .fixtureSampleLogs()
            samplerMock <- new("PopulationSamplerMock", chainSampler=chainSampler )
            sampleLogsOrig <- getSampleLogs(samplerMock)
            sampleLogs(samplerMock) <- fs$sampleLogs
            #getSampleLogs(samplerMock)
            nSampleBefore <- getNSamplePopulations(fs$sampleLogs)
            nSampleBatchPopulations <- c(2L,3L)
            updatedSampler <- setupAndSample(samplerMock
                    , nSample=nSampleBatchPopulations
                    , thin=getThin(chainSampler)
            )
            #
            sampleLogs <- getSampleLogs(updatedSampler)
            parms <- getParametersForPopulation(sampleLogs,2L)
            expect_that( dim(parms), equals(c(length(fs$x0),(nSampleBefore+nSampleBatchPopulations)[2L],getNChainInPopulation(sampleLogs))))
            expect_true( all(is.finite(parms)) )
            pAcceptPops <- getAcceptanceLog(getAcceptanceTracker(updatedSampler))
            expect_that( dim(pAcceptPops), equals(c(getNBlock(sampleLogs),max(nSampleBatchPopulations),getNChain(sampleLogs))))
            expect_true( all(is.finite(pAcceptPops[,min(nSampleBatchPopulations),])) )
            expect_true( all(is.na(pAcceptPops[,-(1:min(nSampleBatchPopulations)),1:2])) )
        })

test_that("sampleLimiting",{
            fs <- .fixtureSampleLogs()
            samplerMock <- new("PopulationSamplerMock", chainSampler=chainSampler )
            sampleLogs(samplerMock) <- fs$sampleLogs
            nSampleBefore <- getNSamplePopulations(fs$sampleLogs)
            updatedSampler <- sampleLimiting(samplerMock)
            #
            sampleLogs <- getSampleLogs(updatedSampler)
            parms <- getParametersForPopulation(sampleLogs,2L)
            nRepresentativeRows  <- computeNRepresentativeRows(getNParameterWithProposal(sampleLogs), getNChainInPopulation(sampleLogs))
            expect_that( dim(parms), equals(c(length(fs$x0),(5L*nRepresentativeRows),getNChainInPopulation(sampleLogs))))
            expect_true( all(is.finite(parms)) )
        })


test_that("newPopulationSampler",{
            sampler <- newPopulationSampler(blocks, x0, diag(0.1*x0)) 
            nSampleBatchPopulations <- c(2L,3L)
            updatedSampler <- setupAndSample(sampler
                    , nSample=nSampleBatchPopulations
                    , thin=2L
            )
            sampleLogs <- getSampleLogs(updatedSampler)
            parms <- getParametersForPopulation(sampleLogs,1L)
            expect_that( dim(parms), equals(c(length(x0),(1L+nSampleBatchPopulations)[1L],getNChainInPopulation(sampleLogs))))
            # here the mock updater is used and initial x0 is not modified
            expect_true( all(parms[,,1] ==x0) )
            for( iChain in 1:getNChainInPopulation(sampleLogs) )
                expect_true( all(parms[,,iChain] == parms[,1,iChain]) )
            updatedSampler2 <- setupAndSample(updatedSampler
                    , nSample=nSampleBatchPopulations
                    , thin=2L
            )
            sampleLogs <- getSampleLogs(updatedSampler2)
            parms <- getParametersForPopulation(sampleLogs,1L)
            expect_that( dim(parms), equals(c(length(x0),(1L+2*nSampleBatchPopulations)[1L],getNChainInPopulation(sampleLogs))))
            # here the mock updater is used and initial x0 is not modified
            expect_true( all(parms[,,1] ==x0) )
        })



