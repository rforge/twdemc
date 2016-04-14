#------------------ Interface for to access Blocks within a chainState
setClass("SampleLogChain", contains="VIRTUAL")

if(!exists("getParameters")) setGeneric("getParameters", function(object) standardGeneric("getParameters"))
if(!exists("getLogDensityComponents")) setGeneric("getLogDensityComponents", function(object) standardGeneric("getLogDensityComponents"))
if(!exists("getProportionAcceptedInInterval")) setGeneric("getProportionAcceptedInInterval", function(object) standardGeneric("getProportionAcceptedInInterval"))

if(!exists("getParametersForSampleIndex")) setGeneric("getParametersForSampleIndex", function(object,iSample) standardGeneric("getParametersForSampleIndex"))
if(!exists("getLogDensityComponentsForSampleIndex")) setGeneric("getLogDensityComponentsForSampleIndex", function(object, iSample) standardGeneric("getLogDensityComponentsForSampleIndex"))
if(!exists("getProportionAcceptedInIntervalForSampleIndex")) setGeneric("getProportionAcceptedInIntervalForSampleIndex", function(object, iSample) standardGeneric("getProportionAcceptedInIntervalForSampleIndex"))
