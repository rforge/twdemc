setClass("JumpProposer", contains="VIRTUAL" )

if(!exists("adjustToAcceptanceRate")) setGeneric("adjustToAcceptanceRate", 
            function(object, acceptanceTracker,...) 
                standardGeneric("adjustToAcceptanceRate"))

if(!exists("proposeJumps")) setGeneric("proposeJumps", 
            function(object, nGeneration, sampleLogs, iPopulationsStep, iCurrentSamples,...) 
                standardGeneric("proposeJumps"))

if(!exists("initializeDimensionalSettings")) setGeneric("initializeDimensionalSettings", 
            function(object,...) 
                standardGeneric("initializeDimensionalSettings"))

# proposals may be tailored during burnin
if(!exists("isBurnin<-")) setGeneric("isBurnin<-", function(object,value) standardGeneric("isBurnin<-"))
