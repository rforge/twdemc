#' @include SettingsWithDefaults.R

#' @export
setClass("AcceptanceTrackerSettings", contains="SettingsWithDefaults")

#' @export
setMethod("initializeDefaults", signature="AcceptanceTrackerSettings", 
        function(object
            ### Set default settings.
        ){
            ##value<< Settings object with the following entries set:
            object@ctrl <- within(object@ctrl,{        
                        widthAcceptWindow = 5L          ##<< scalar integer: how many intervals are to be averaged to estimate acceptance rate  
                        initialAcceptanceRate = 0.25    ##<< acceptance rate assumed initially, when no acceptance rates have been recorded yet
                        ##end<<
                    })
            return(object)
        })


