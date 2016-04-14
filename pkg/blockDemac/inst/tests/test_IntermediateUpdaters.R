#require(testthat)
context("IntermediateUpdaters")

theta0 <- c(a=2,b=3)

intermediateSpecs <- list(
        model1 = intermediateSpec(
                function(theta, intermediates){ pred <- theta*10 }
                ,parameters = names(theta0)
        )
        ,int2 = intermediateSpec(
                function(theta, intermediates){ intermediates[["model1"]]*2 }
                ,intermediates = "model1"
        )
)

logDensityComponents0 =c(m1=-100, m2=-100, parms=-200)
#intermediates0 = lapply( intermediateSpecs, function(entry){list()}) 
chainState <- new("ChainState", parameters=theta0, logDensityComponents=logDensityComponents0
        ,intermediateIds = names(intermediateSpecs))
#,intermediates = intermediates0)


intUpdaters <- newIntermediateUpdaters(intermediateSpecs, names(theta0))

test_that("initializeIntermediateUpdaters",{
            intUpdater <- getIntermediateUpdater(intUpdaters,"model1")
            expect_that( getParametersUsed(intUpdater), equals(names(theta0)))
            expect_that( getBlockIndicesIParametersUsed(intUpdater), equals(seq_along(theta0)))
            expect_that( getBlockIndicesIntermediateIdsUsed(intUpdater), equals( character(0)))
            intUpdater <- getIntermediateUpdater(intUpdaters,"int2")
            expect_that( getParametersUsed(intUpdater), equals(character(0)))
            expect_that( getBlockIndicesIParametersUsed(intUpdater), equals(integer(0)))
            expect_that( getBlockIndicesIntermediateIdsUsed(intUpdater), equals("model1"))            
            expect_that( intUpdater@requisiteUpdaterIds, equals("model1"))            
        })

test_that("initialize with unknown parameters",{
            intermediateSpecs <- list(
                    model1 = intermediateSpec(
                            function(theta, intermediates){ pred <- theta*10 }
                            ,parameters = "nonExistingParameter"
                    )
            )
            expect_error(
                    intUpdaters <- newIntermediateUpdaters(intermediateSpecs, names(theta0))
            )
        })

test_that("error on non-unique intermediateSpecification names",{
            intSpecs <- list(
                    B= blockSpecMock(c("a"))
                    ,B = blockSpecMock(c("b"))
            )
            expect_error(
                    intUpdaters <- newIntermediateUpdaters(intSpecs, c("a","b"))
            )
            intSpecs <- list(
                    blockSpecMock(c("a"))
                    ,B = blockSpecMock(c("b"))
            )
            expect_error(
                    intUpdaters <- newIntermediateUpdaters(intSpecs, c("a","b"))
            )
            intSpecs <- list(
                    blockSpecMock(c("a"))
            )
            expect_error(
                    intUpdaters <- newIntermediateUpdaters(intSpecs, c("a","b"))
            )
        })

test_that("update and get Intermediate in ChainState",{
            intUpdater <- getIntermediateUpdater(intUpdaters,"model1")
            intermediate <- getIntermediateFromUpdater( intUpdater, chainState, isWarnOnEmpty=FALSE )
            expect_that( intermediate, equals(list()))            
            updatedChainState <- letUpdaterUpdateIntermediateInChainState(intUpdater, chainState, intUpdaters@intermediateUpdaters)            
            expect_that( getIntermediateFromUpdater( intUpdater, updatedChainState ), equals(theta0*10))
            #
            intUpdater <- getIntermediateUpdater(intUpdaters,"int2")
            intermediate <- getIntermediateFromUpdater( intUpdater, chainState, isWarnOnEmpty=FALSE )
            expect_that( intermediate, equals(list()))            
            updatedChainState <- letUpdaterUpdateIntermediateInChainState(intUpdater, chainState, intUpdaters@intermediateUpdaters)            
            expect_that( getIntermediateFromUpdater( intUpdater, updatedChainState ), equals(theta0*10*2))
            # check that the intermediate which this state depends on also has been updated 
            expect_that( getIntermediateFromUpdater( getIntermediateUpdater(intUpdaters,"model1"), updatedChainState ), equals(theta0*10))
        })

test_that("get updaters depending on intermediateIds",{
            depIntUpdaters <- getIntermediateUpdatersDependingOnIntermediates(intUpdaters, character(0))
            expect_that( length(depIntUpdaters), equals(0L))
            depIntUpdaters <- getIntermediateUpdatersDependingOnIntermediates(intUpdaters,"int2")
            expect_that( length(depIntUpdaters), equals(0L))
            depIntUpdaters <- getIntermediateUpdatersDependingOnIntermediates(intUpdaters,"model1")
            expect_that( names(depIntUpdaters), equals("int2"))
        })


test_that("get updaters depending on parameters",{
            depIntUpdaters <- getDependentIntermediateUpdaters(intUpdaters,5L)
            expect_that( length(depIntUpdaters), equals(0L))
            depIntUpdaters <- getDependentIntermediateUpdaters(intUpdaters,1L)
            expect_that( sort(names(depIntUpdaters)), equals(sort(c("model1","int2"))) )
        })





