#require(testthat)
context("BlockUpdaters")


fUpdateBlockAddTerm <- function(
        ### update function that appends a logString to variable log in environment env
        x, iParametersToUpdate
        , lowerParBounds, upperParBounds
        , intermediates
        , term
        , env
        , logString="fUpdaterBlockAddTerm"
){
    if( !missing(env) ) env$log <- paste0(env$log, logString)
    list(
            isUpdated = TRUE 
            ,xC = x[ iParametersToUpdate ] + term
            ,intermediate = NULL    # actually dont return NULL, better list(), here testing handled
    )
}

x0 <- c(a=1, b=2, c=3)

blocks <- list(
        add7 = blockSpec("a",names(x0), new("FunctionBasedBlockUpdater", 
                        fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=7))), 
        add9 = blockSpec(c("b","c"),names(x0),new("FunctionBasedBlockUpdater",
                        fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=9)))
)

blocksSingle <- list(
        all = blockSpec(,names(x0), new("FunctionBasedBlockUpdater", 
                        fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=7)))
)

blocksMetr <- list(
        lden1 = blockSpec("a",names(x0), new("MetropolisBlockUpdater", 
                        fLogDensity=function(x,intermediates,logDen){c(obs1=-1)},
                        logDensityComponentNames = c("obs1"))), 
        lden23 = blockSpec(c("b","c"),c("b","c"),new("MetropolisBlockUpdater", 
                        fLogDensity=function(x,intermediates,logDen){c(obs2=-2,parms=-3)},
                        logDensityComponentNames = c("obs2","parms")))
)

intermediateSpecs <- list()


logEnv <- new.env(parent=emptyenv())
logEnv$log <- ""

# chain of dependent intermediates that add a logString to env$log
# (a,b) -> intab -> intintab
# (c) -> intc
intermediateSpecsDep <- list(
        intab = intermediateSpec(
                function(theta, intermediates, env){
                    if( !missing(env) ) env$log <- paste0(env$log,",intab")
                    c(intab=theta["a"]+theta["b"]) 
                }
                ,argsUpdateFunction=list(env=logEnv)
                ,parameters = c("a","b")
        )
        ,intintab = intermediateSpec(
                function(theta, intermediates, env){
                    if( !missing(env) ) env$log <- paste0(env$log,",intintab")
                    c(intintab=intermediates[[1]]*2) 
                }
                ,argsUpdateFunction=list(env=logEnv)
                ,intermediates = "intab"
        )
        ,intc = intermediateSpec(
                function(theta, intermediates, env){
                    if( !missing(env) ) env$log <- paste0(env$log,",intc")
                    c(intc=theta*2) 
                }
                ,argsUpdateFunction=list(env=logEnv)
                ,parameters = c("c")
        )
)

# blocks dependent on above intermediates
# see dependencyGraph in inst/develop/test_BlockUpdatersDep.graphml
# (a) -> A -> a
# (intab) -> B -> b
# (intintab,intc) -> C -> c
blocksDep <- list(
        A = blockSpec("a",, new("FunctionBasedBlockUpdater", 
                        fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=7,env=logEnv,logString=",A"))) 
        ,B = blockSpec("b",, new("FunctionBasedBlockUpdater",
                        fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=9,env=logEnv,logString=",B"))
                ,intermediatesUsed=c("intab") )
        ,C = blockSpec("c",, new("FunctionBasedBlockUpdater",
                        fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=9,env=logEnv,logString=",C"))
                ,intermediatesUsed=c("intintab","intc") )
)

test_that("access blockSpecifications",{
            # testing correct access of BlockSpecifications with direct dependencies
            expect_that( getParametersUsed(blocks$add7), equals(names(x0)) )
            #
            expect_that( getParametersUsed(blocksDep$A), equals("a") )
            expect_that( getIntermediatesUsed(blocksDep$A), equals(character(0)) )
            expect_that( getParametersUsed(blocksDep$B), equals("b") )
            expect_that( getIntermediatesUsed(blocksDep$B), equals(c("intab")) )
            expect_that( getParametersUsed(blocksDep$C), equals("c") )
            expect_that( getIntermediatesUsed(blocksDep$C), equals(c("intintab","intc")) )
        })

test_that("initialization",{
            blockUpdaters <- newBlockUpdaters(blocks, intermediateSpecs, names(x0))
            expect_that( getLength(blockUpdaters), equals(2) )
            expect_that( getNames(blockUpdaters), equals(c("add7","add9")) )
            blockUpdaters <- newBlockUpdaters(blocksMetr, intermediateSpecs, names(x0))
            blockDimensions <- getBlockDimensions(blockUpdaters)        
            expect_that( getLogDensityComponentNames(blockDimensions), equals(c("obs1","obs2","parms")) )        
        })


test_that("error on initialization nBlock=0",{
            expect_error(
                    blockUpdaters <- newBlockUpdaters(list(), intermediateSpecs, names(x0))
            )
        })

test_that("error on non-unique blockSpecification names",{
            blocks <- list(
                    B= blockSpecMock(c("a"))
                    ,B = blockSpecMock(c("b"))
            )
            expect_error(
                    blockUpdaters <- newBlockUpdaters(blocks, list(), c("a","b"))
            )
            blocks <- list(
                    blockSpecMock(c("a"))
                    ,B = blockSpecMock(c("b"))
            )
            expect_error(
                    blockUpdaters <- newBlockUpdaters(blocks, list(), c("a","b"))
            )
            blocks <- list(
                    blockSpecMock(c("a"))
            )
            expect_error(
                    blockUpdaters <- newBlockUpdaters(blocks, list(), c("a","b"))
            )
        })


test_that("error on parameters double updated",{
            blocks <- list(
                    A = blockSpecMock("a")
                    ,B = blockSpecMock(c("a","b"))
            )
            expect_error(
                    blockUpdaters <- newBlockUpdaters(blocks, list(), c("a","b"))
            )
            blocks <- list(
                    A = blockSpecMock() # empty defaults to all parameters updated
                    ,B = blockSpecMock(c("b"))
            )
            expect_error(
                    blockUpdaters <- newBlockUpdaters(blocks, list(), c("a","b"))
            )
        })

test_that("error on updating non-existing parameters",{
            blocks <- list(
                    B = blockSpecMock(c("a","b"))
            )
            expect_error(
                    blockUpdaters <- newBlockUpdaters(blocks, list(), c("a"))
            )
        })

test_that("error on parameters not updated by block",{
            blocks <- list(
                    B = blockSpecMock(c("a"))
            )
            expect_error(
                    blockUpdaters <- newBlockUpdaters(blocks, list(), c("a","b"))
            )
        })



test_that(".getRequisiteIntermediateUpdaters",{
            # check that recursively got all intermediateUpdater necessary fro blockUpdater C
            blockUpdaters <- newBlockUpdaters(blocksDep, intermediateSpecsDep, names(x0))
            bUpdater <- getBlockUpdaterForName(blockUpdaters,"C")
            intUpdaters <- .getRequisiteIntermediateUpdaters(blockUpdaters, bUpdater)
            expect_that( length(intUpdaters), equals(3L) )
            expect_that( sort(names(intUpdaters)), equals(c("intab", "intc", "intintab")))
        })

test_that(".getDependentIntermediateUpdaters",{
            # check that recursively got all IntermedidateUpdater depending on parameters (b) updated by Block B
            blockUpdaters <- newBlockUpdaters(blocksDep, intermediateSpecsDep, names(x0))
            bUpdater <- getBlockUpdaterForName(blockUpdaters,"B")
            intUpdaters <- .getDependentIntermediateUpdaters(blockUpdaters, bUpdater)
            expect_that( length(intUpdaters), equals(2L) )
            expect_that( sort(names(intUpdaters)), equals(c("intab", "intintab")))
        })

test_that("getDependentBlockUpdaters",{
            # check that recursively got all Blockupdaters depending on parameters (a) updated by Block A
            blockUpdaters <- newBlockUpdaters(blocksDep, intermediateSpecsDep, names(x0))
            bUpdater <- getBlockUpdaterForName(blockUpdaters,"A")
            depUpdaters <- getDependentBlockUpdaters(blockUpdaters,"A")
            expect_that( length(depUpdaters), equals(3L) )
            parametersToUpdate <- sapply(depUpdaters, function(depUpdater){ depUpdater@iParametersToUpdate })
            expect_that( sort(names(parametersToUpdate)), equals(c("a","b", "c")))
        })



test_that("intermediateUpdatersInitialization",{
            # see dependencyGraph in inst/develop/test_BlockUpdatersDep.graphml
            blockUpdaters <- newBlockUpdaters(blocksDep, intermediateSpecsDep, names(x0))
            intermediateUpdaters <- getIntermediateUpdaters(blockUpdaters)
            expect_that( getNames(intermediateUpdaters), equals(names(intermediateSpecsDep)))
            #
            # track dependencies in blockUpdaters
            expect_that( blockUpdaters@dependentIntermediateIds$A, equals(c("intab", "intintab")) )
            expect_that( blockUpdaters@dependentIntermediateIds$B, equals(c("intab", "intintab")) )
            expect_that( blockUpdaters@dependentIntermediateIds$C, equals(c("intc")) )
            expect_that( blockUpdaters@dependentBlockIndices$A, equals(1:3) )
            expect_that( blockUpdaters@dependentBlockIndices$B, equals(2:3) )
            expect_that( blockUpdaters@dependentBlockIndices$C, equals(3L) )
        })

test_that("computeUpdatedChainState",{
            blockUpdaters <- newBlockUpdaters(blocksDep, intermediateSpecsDep, names(x0))
            subSpace(blockUpdaters) <- initializeIndexedBoundsForParameterNames(new("SubSpace"), names(x0))
            chainState <- new("ChainState" 
                    ,parameters=x0 
                    ,logDensityComponents=c(obs1=-100, obs2=-200)
                    ,intermediateIds=names(intermediateSpecsDep)
            )
            stepInfo <- new("StepInfoImpl", step = c(a=1,b=2,c=3), nLogDensityComponents=2 )
            logEnv$log <- ""            
            updatedChainState <- computeUpdatedChainState(chainState, blockUpdaters, stepInfo)
            expect_that(logEnv$log, equals(",A,intab,B,intab,intintab,intc,C") )
            logEnv$log <- ""            
            updatedChainState2 <- computeUpdatedChainState(updatedChainState, blockUpdaters, stepInfo)
            #logEnv$log
            expect_that(logEnv$log, equals(",A,intab,B,intab,intintab,intc,C") )
        })

test_that("computeUpdatedChainState modelSigma",{
            intermediateSpecsDep2 <- list(
                    model1 = intermediateSpec(
                            function(theta, intermediates, env){
                                if( !missing(env) ) env$log <- paste0(env$log,",model1")
                                c(intab=theta["a"]+theta["b"]) 
                            }
                            ,argsUpdateFunction=list(env=logEnv)
                            ,parameters = c("a","b")
                    )
            )
            blocksDep2 <- list(
                    obs = blockSpec(c("a","b"),, new("MetropolisBlockUpdater", 
                                    fLogDensity=function(x,intermediates,logDen, env){
                                        if( !missing(env) ) env$log <- paste0(env$log,",obs")
                                        c(obs1=-1)
                                    }
                                    ,argsFLogDensity=list(env=logEnv)
                                    ,logDensityComponentNames = c("obs1")) 
                            ,intermediatesUsed=c("model1") )
                    ,sigma = blockSpec("c",, new("FunctionBasedBlockUpdater",
                                    fUpdateBlock=fUpdateBlockAddTerm, argsFUpdateBlock=list(term=9,env=logEnv,logString=",sigma"))
                            ,intermediatesUsed=c("model1") )
            )
            blockUpdaters <- newBlockUpdaters(blocksDep2, intermediateSpecsDep2, names(x0))
            subSpace(blockUpdaters) <- initializeIndexedBoundsForParameterNames(new("SubSpace"), names(x0))
            chainState <- new("ChainState" 
                    ,parameters=x0 
                    ,logDensityComponents=c(obs1=-100)
                    ,intermediateIds=names(intermediateSpecsDep2)
            )
            stepInfo <- new("StepInfoImpl", step = c(a=1,b=2,c=3), nLogDensityComponents=2 )
            logEnv$log <- ""            
            updatedChainState <- computeUpdatedChainState(chainState, blockUpdaters, stepInfo)
            # model1 initial state, model1 proposed state, obs, sigma
            expect_that(logEnv$log, equals(",model1,model1,obs,sigma") )
            logEnv$log <- ""            
            updatedChainState2 <- computeUpdatedChainState(updatedChainState, blockUpdaters, stepInfo)
            # model1 proposed state, obs, sigma
            expect_that(logEnv$log, equals(",model1,obs,sigma") )
        })




test_that("getBlockUpdaterFor on existing entry",{
            bl <- newBlockUpdaters(blocks, intermediateSpecs, names(x0))
            blockUpdater <- getBlockUpdaterForName(bl,"add7")
            expect_that(blockUpdater@argsFUpdateBlock, equals(list(term=7)))
            blockUpdater <- getBlockUpdaterForIndex(bl,2L)
            expect_that(blockUpdater@argsFUpdateBlock, equals(list(term=9)))
        })

test_that("getBlockUpdaterFor on several or no entries yield error",{
            bl <- newBlockUpdaters(blocks, intermediateSpecs, names(x0))
            expect_error(
                    blockUpdater <- getBlockUpdaterForName(bl,c("add7","add9"))
            )
            expect_error(
                    blockUpdater <- getBlockUpdaterForName(bl,character(0))
            )
        })

test_that("getBlockUpdaterFor on nonexisting entry yield error",{
            bl <- newBlockUpdaters(blocks, intermediateSpecs, names(x0))
            expect_error(
                    blockUpdater <- getBlockUpdaterForName(bl,c("nonExisting"))
            )
        })

test_that("initializeBlockRelations",{
            bl <- newBlockUpdaters(blocks, intermediateSpecs, names(x0))
            #bl <- initializeBlockRelations(bl, blocks, names(x0) )
            expect_that( getIParametersThatRequireJumpProposal(bl), equals(integer(0)))
            add7 <- getBlockUpdaterForName(bl,"add7")
            expect_that( getBlockIndicesIParametersToUpdate(add7), is_equivalent_to(1))
            expect_that( getBlockIndicesIParametersUsed(add7), is_equivalent_to(1:3))
            expect_that( getBlockIndicesILogDensityComponents(add7), is_equivalent_to(integer(0)))
            #
            bl <- newBlockUpdaters(blocksSingle, intermediateSpecs, names(x0))        
            #bl <- initializeBlockRelations(bl, blocksSingle, names(x0) )
            expect_that( getIParametersThatRequireJumpProposal(bl), equals(integer(0)))
            blockUpdater <- getBlockUpdaterForName(bl,"all")
            expect_that( getBlockIndicesIParametersToUpdate(blockUpdater), is_equivalent_to(1:3))
            expect_that( getBlockIndicesIParametersUsed(blockUpdater), is_equivalent_to(1:3))
            expect_that( getBlockIndicesILogDensityComponents(blockUpdater), is_equivalent_to(integer(0)))
            #
            bl <- newBlockUpdaters(blocksMetr, intermediateSpecs, names(x0))
            #bl <- initializeBlockRelations(bl, blocksMetr, names(x0) )
            expect_that( getIParametersThatRequireJumpProposal(bl), equals(structure(1:3, .Names = c("a", "b", "c"))))
            blockUpdater <- getBlockUpdaterForName(bl,"lden1")
            expect_that( getBlockIndicesIParametersToUpdate(blockUpdater), is_equivalent_to(1))
            expect_that( getBlockIndicesIParametersUsed(blockUpdater), is_equivalent_to(1:3))
            expect_that( getBlockIndicesILogDensityComponents(blockUpdater), is_equivalent_to(1L))
            blockUpdater <- getBlockUpdaterForName(bl,"lden23")
            expect_that( getBlockIndicesIParametersToUpdate(blockUpdater), is_equivalent_to(2:3))
            expect_that( getBlockIndicesIParametersUsed(blockUpdater), is_equivalent_to(2:3))
            expect_that( getBlockIndicesILogDensityComponents(blockUpdater), is_equivalent_to(2:3))
        })

test_that("getBlockDimensions",{
            bl <- newBlockUpdaters(blocks, intermediateSpecs, names(x0))
            blockDimensions <- getBlockDimensions(bl)
            expect_that( getParameterNames(blockDimensions), equals(names(x0)) )
        })












