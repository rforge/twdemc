setClass("SampleDimensionsImpl",
        representation(
                nSamplePopulations = "integer"
                ,nChainInPopulation = "integer"
        )
        ,prototype(nSamplePopulations=4L, nChainInPopulation=3L)
)

setGeneric("getIPopulationsForChains", function(object, iChain) standardGeneric("getIPopulationsForChains"))
setMethod("getIPopulationsForChains", signature=c("SampleDimensionsImpl","integer"), 
        function(object, iChain) {
            ((iChain-1L) %/% object@nChainInPopulation)+1L 
        })


s1 <- new("SampleDimensionsImpl")
getIPopulationsForChains(s1,4L)




