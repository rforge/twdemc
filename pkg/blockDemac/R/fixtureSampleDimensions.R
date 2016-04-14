.fixtureSampleDimensions <- function( 
        nSample = c(26L,30L)      # acommodate 20 states of pre-History
        ,nChainInPopulation = 2L
        ,nChainGroupsInPopulation = 1L #max(1L,as.integer(nChainInPopulation %/% 3L))
        ,countAcceptedInInterval0 = c(met1=2,met2=3) # number of accepted steps at first sample, ie. thinning interval
        ,blockSpecifications = list(
                m1=blockSpecMock(c("a"),,c("obs1"))
                ,m2Parms =blockSpecMock(c("b","c"),,c("obs2","parms"))
                )
){
    parNames <- do.call(c, structure(lapply(blockSpecifications, getParametersToUpdate), names=NULL))
    blockDimensions <- newBlockDimensions(blockSpecifications, parNames )
    sampleDimensions <- initializeSampleDimensionsImpl( new("SampleDimensionsImpl")
        ,blockDimensions
        ,nSamplePopulations = nSample
        ,nChainInPopulation = nChainInPopulation
        ,nChainGroupsInPopulation = nChainGroupsInPopulation
    )
    attr(sampleDimensions,"blockSpecifications") <- blockSpecifications
    sampleDimensions
}
attr(.fixtureSampleDimensions,"ex") <- function(){
    .fixtureSampleDimensions()
}
