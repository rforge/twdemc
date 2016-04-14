.computeExponentialDecreaseVector <- function(
        ### Calculates the temperature for an exponential decrease from \code{TStart} to \code{Tend} after \code{nGen} steps. 	
        TStart			##<< numeric vector: the initial temperature (before the first step at iGen=0)
        , TEnd=1	##<< numeric vector: the temperature at the last step
        , nGen		##<< integer scalar: the number of genrations	
){
    mapply(.computeExponentialDecreaseScalar, TStart=TStart, TEnd=TEnd, nGen=nGen)
}
attr(.computeExponentialDecreaseVector,"ex") <- function(){
    res <- .computeExponentialDecreaseVector(TStart=c(20,10),TEnd=c(2,1),nGen=100) 
    matplot( 1:100, res, type="l" )	
    res <- .computeExponentialDecreaseVector(TStart=c(20,10),TEnd=c(2),nGen=100) 
    matplot( 1:100, res, type="l" )	
}


.computeExponentialDecreaseScalar <- function( 
        ### Calculates the temperature for an exponential decrease from \code{TStart} to \code{Tend} after \code{nGen} steps. 	
        TStart			##<< numeric scalar: the initial temperature (before the first step at iGen=0)
        , TEnd=1	##<< numeric scalar: the temperature at the last step
        , nGen		##<< integer scalar: the number of genrations	
){
    # calcDEMCTemp
    ##seealso<< 
    ## \code{\link{calcTemperatedLogDen.matrix}}
    ## \code{\link{twDEMC}}
    if( nGen < 1) return( numeric(0) )
    b = TStart
    a = log(TEnd/TStart)/nGen
    iGen=1:nGen ##<< integer vector: the steps for which to calculate the Temperature	
    ##value<< vector of Temperatures corresponding to steps iGen
    b*exp( a*iGen )
}
attr(.computeExponentialDecreaseScalar,"ex") <- function(){
    plot( 1:100, .computeExponentialDecreaseScalar(TStart=100,TEnd=5,nGen=100) )	
}
