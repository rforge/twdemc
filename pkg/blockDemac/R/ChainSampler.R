setClass("ChainSampler", contains="VIRTUAL" )

if(!exists("getBlockDimensions")) setGeneric("getBlockDimensions", function(object) standardGeneric("getBlockDimensions"))
if(!exists("getIParametersThatRequireJumpProposal")) setGeneric("getIParametersThatRequireJumpProposal", function(object) standardGeneric("getIParametersThatRequireJumpProposal"))

if(!exists("setRangeSpecs")) setGeneric("setRangeSpecs", function(object, ...) standardGeneric("setRangeSpecs"))
if(!exists("subSpace<-")) setGeneric("subSpace<-", function(object,value) standardGeneric("subSpace<-"))
if(!exists("thin<-")) setGeneric("thin<-", function(object,value) standardGeneric("thin<-"))

if(!exists("getBlockIndicesIntermediateIdsUsed")) setGeneric("getBlockIndicesIntermediateIdsUsed", function(object,...) standardGeneric("getBlockIndicesIntermediateIdsUsed"))
if(!exists("computeChainState")) setGeneric("computeChainState", function(object, ...) standardGeneric("computeChainState"))
if(!exists("sampleRange")) setGeneric("sampleRange", function(object) standardGeneric("sampleRange"))

if(!exists("getChainState")) setGeneric("getChainState", function(object) standardGeneric("getChainState"))
if(!exists("getSampleLog")) setGeneric("getSampleLog", function(object) standardGeneric("getSampleLog"))

