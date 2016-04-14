#' @include BlockDimensions.R

#' @export
setClass("SampleDimensions", contains="BlockDimensions")

if(!exists("getNSamplePopulations")) setGeneric("getNSamplePopulations", function(object) standardGeneric("getNSamplePopulations"))
if(!exists("getNPopulation")) setGeneric("getNPopulation", function(object) standardGeneric("getNPopulation"))
if(!exists("getNChain")) setGeneric("getNChain", function(object) standardGeneric("getNChain"))
if(!exists("getNChainInPopulation")) setGeneric("getNChainInPopulation", function(object) standardGeneric("getNChainInPopulation"))
if(!exists("getNChainGroupsInPopulation")) setGeneric("getNChainGroupsInPopulation", function(object) standardGeneric("getNChainGroupsInPopulation"))
if(!exists("getIChainsForPopulation")) setGeneric("getIChainsForPopulation", function(object, iPopulation) standardGeneric("getIChainsForPopulation"))
#if(!exists("getIPopulationChains")) setGeneric("getIPopulationChains", function(object) standardGeneric("getIPopulationChains"))
if(!exists("getIChainsForPopulations")) setGeneric("getIChainsForPopulations", function(object, iPopulations) standardGeneric("getIChainsForPopulations"))
if(!exists("getIPopulationsForChains")) setGeneric("getIPopulationsForChains", function(object, iChain) standardGeneric("getIPopulationsForChains"))




