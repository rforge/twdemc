# slowly decreasing Temperature (cost reduction factor)

R.methodsS3::setMethodS3("calcTemperatedLogDen","default", function(
                ### Rescale Log-Densities by given Temperatures
                x			##<< numeric matrix (nStep x nResComp): logDensities for each component at each step
                ,temp		##<< numeric vector (nResComp): Temperature, i.e. cost reduction factor
                ,...
        ){
            # calcTemperatedLogDen.default
            ##seealso<<   
            ## \code{\link{twDEMC}}
            
            ##details<< 
            ## There are several function to help with temperating log-Densities
            ## \itemize{
            ## \item{ Rescale Log-Densities by given Temperatures: this function  } 
            ## \item{ Rescale Log-Densities within chains \code{\link{calcTemperatedLogDenChains.array}} }
            ## \item{ Generating exponential temperature series: \code{\link{calcDEMCTemp}}  }
            ## \item{ Calculating a base temperature from stream temperatures: \code{\link{calcBaseTemp}}  }
            ## \item{ Calculating a base temperature from cost values: \code{\link{calcBaseTempSk}}  }
            ## \item{ Calculating a Stream temperatures from base temperature: \code{\link{calcStreamTemp}}  }
            ##}
            # # \item{ Rescale Log-Densities by given Temperatures with an offst: \code{\link{calcTemperatedLogDenOffset.default}}  }
            
            if( ncol(x) == 1){
                x/temp
            }else{
                t(t(x)/temp)        # transform necessary to match nResComp of x to temp 
            }
            ##value<< numeric matrix (nStep x nResComp), rescaled logDen
            
        })
attr(calcTemperatedLogDen,"ex") <- function(){
    data(twdemcEx1)
    logDen0 <- stackChains(concatPops(twdemcEx1)$resLogDen)
    refLogDen <- apply(logDen0,2,max)
    temp = c(obs=20, parms=1)
    logDen <- t( (t(logDen0)-refLogDen)*temp +refLogDen )
    
    logDenT <- calcTemperatedLogDen(logDen,temp)
    .exp <- logDen0
    names(dimnames(.exp)) <- NULL 
    names(dimnames(logDenT)) <- NULL 
    all.equal( .exp, logDenT )
}



calcBaseTemp <- function(
        ### calculate Base temperature from given temperatures of data streams and corresponding number of observations.
        temp        ##<< numeric vector (nResultComp): temperature, i.e. variance inflation factor
        , nObs      ##<< integer vector (nResultComp): number of observations
        , TFix=rep( NA_real_, length(temp))       ##<< numeric vector (nResultComp): fixed temperature for components, non-finite for those with varying temperature
        ## , alternatively a named vector listing only the components with fixed temperatures (temp must have names then too)
        , iNonFixTemp=which(!is.finite(TFix))     ##<< integer vector: index of result components, which Temperature is fixed 
){
    ##seealso<< \link{calcStreamTemp}, \link{twDEMC}
    ##seealso<< \code{\link{calcTemperatedLogDen.default}}
    
    ##details<<
    ## The Temperature, i.e. variance inflation factor, scales with the number of observations by
    ## T = 1 + (T0 - 1) * n.  
    ## A common base temperature is the maximium across calculated stream base temperatures
    if( length(nObs) != length(temp) )
        nObs <- nObs[names(temp)]
    if( length(TFix) != length(temp) ){
        TFixN <- structure( rep( NA_real_, length(temp)), names=names(temp) )
        TFixN[names(TFix)] <- TFix
        TFix <- TFixN
    }
    tempBaseK <- (temp -1)/nObs + 1    # one estimate for base factor
    ##value<< numeric scalar: The temperature, i.e. variance inflation factor, for a single observation.
    tempBase <- max( tempBaseK[iNonFixTemp], na.rm=TRUE)    # take the maximum of the bases, for components that do not have fixed temperature
    tempBase
}
attr(calcBaseTemp,"ex") <- function(){
    data(twdemcEx1)   
    .nObs <- c(parms=getNParms(twdemcEx1), obs=length(twdemcEx1$dInfos[[1]]$argsFLogDen$obs) )
    .T <- getCurrentTemp(twdemcEx1)
    calcBaseTemp( .T, .nObs[names(.T)], TFix=c(parms=1) )
}

#trace(calcStreamTemp, recover)  #untrace(calcStreamTemp)
calcStreamTemp <- function(
        ### scale base temperature to given number of observations.
        tempBase    ##<< numeric scalar: The temperature, i.e. variance inflation factor, for a single observation.
        , nObs      ##<< integer vector (nResultComp): number of observations
        , TFix=rep( NA_real_, length(nObs))      ##<< numeric vector (nResultComp): fixed temperature for components, non-finite for those with varying temperature
        , iFixTemp=which(is.finite(TFix))        ##<< integer vector: index of result components, which Temperature is fixed 
){
    ##seealso<< \link{calcBaseTemp}, \link{twDEMC}
    ##seealso<< \code{\link{calcTemperatedLogDen.default}}
    
    if( length(TFix) != length(nObs) ){
        TFixN <- structure( rep( NA_real_, length(nObs)), names=names(nObs) )
        TFixN[names(TFix)] <- TFix
        TFix <- TFixN
    }
    temp <- 1 +(tempBase-1)*nObs
    temp[iFixTemp] <- TFix[iFixTemp]
    ##value<< numeric vector: The Temperatures for the observation streams.
    temp
}

calcBaseTempSk <- function(
        ### Calculate the mean variance factor from given cost and observation numbers.        
        Sk               ##<< numeric matrix (nCases, nResultComp): misfit (=-2*logDensity) per result components
        , nObs=1         ##<< number of observations per result component 
        , TFix=rep( NA_real_, ncol(Sk))        ##<< numeric vector (nResultComp): fixed temperature for components, non-finite for those with varying temperature
        , iFixTemp=which(is.finite(TFix))        ##<< integer vector: index of result components, which Temperature is fixed 
        , iNonFixTemp=which(!is.finite(TFix))    ##<< integer vector: index of result components, which Temperature is fixed
        , isVerbose = FALSE ##<< set to TRUE to report mean base temperatures per stream
){
    ##seealso<< \code{\link{calcTemperatedLogDen.default}}
    if( length(iNonFixTemp) ){
        # aling nObs and 
        Sknf <- Sk[,iNonFixTemp ,drop=FALSE]
        nObsA <- nObs[ colnames(Sknf) ]
        T0sr <- T0s <- t(Sknf)/nObsA
        T0s[] <- pmax(1, T0sr)  # [] in order to keep attributes
        if( isTRUE(isVerbose)) print(signif(rowMeans(T0s)-1,2))
        T0 <- mean(T0s)        
    }else{
        1
    }
    ##value<< numeric scalar: The temperature, i.e. variance inflation factor, for a single observation.
}

R.methodsS3::setMethodS3("calcTemperatedLogDen","twDEMCPops", function(
	### Rescale Log-Densities by given Temperatures
	x							##<< object of class twDEMCPops
	,temp=getCurrentTemp(x)		##<< numeric vector (nResComp): Temperature, i.e. cost reduction factor
	,...
){
    ##seealso<<   \code{\link{calcTemperatedLogDen.default}}
    ##seealso<<   \code{\link{calcTemperatedLogDenChains.twDEMC}}
    logDen <- stackChains(concatPops(x)$resLogDen)
	calcTemperatedLogDen(logDen,temp,...)
    ##value<< array (nCases * nComp) of stacked rescaled logDensities
})

R.methodsS3::setMethodS3("calcTemperatedLogDenChains","array", function(
		### Rescale Log-Densities by given Temperatures within chains
		x					##<< numeric matrix (nStep x nResComp x nChain), e.g. pop$resLogDen
		,temp				##<< numeric vector (nResComp): Temperature, i.e. cost reduction factor
		,...
	){
        ##seealso<<   \code{\link{calcTemperatedLogDen.default}}
        #calcTemperatedLogDenChains.array
		#refLogDen=apply(x,2,max)		# same minimum across all chains
		#refLogDen=apply(x,2,quantile,probs=0.9)		# same reference across all chains
		#iChain <- 1
		rL <- lapply( 1:dim(x)[3], function(iChain){ 
				#calcTemperatedLogDen.default( adrop(x[,,iChain ,drop=FALSE],3), temp, refLogDen=refLogDen )
                calcTemperatedLogDen.default( adrop(x[,,iChain ,drop=FALSE],3), temp )
            })
		logDenT <- abind(rL, rev.along=0)
		##value<< numeric array (nStep x nComp x nChain): rescaled Log-Density for each chain.
		##seealso<< \code{\link{calcTemperatedLogDen.default}}
	})


R.methodsS3::setMethodS3("calcTemperatedLogDenChains","twDEMC", function(
        ### Rescale Log-Densities by given Temperatures within chains
		x							##<< object with entry resLogDen, a numeric matrix (nStep x nResComp x nChain), e.g. twDEMC or pop in twDEMCPops
		,temp=getCurrentTemp(x)		##<< numeric vector (nResComp): Temperature, i.e. cost reduction factor
		,...
	){
        ##seealso<<   \code{\link{calcTemperatedLogDen.default}}
        #calcTemperatedLogDenChains.twDEMC
		calcTemperatedLogDenChains.array( x$resLogDen, temp, ... )
        # array of (nCases * nComp * nChains) of unstacked rescaled logDensity components
	})


calcDEMCTemp <- function( 
	### Calculates the temperature for an exponential decrease from \code{T0} to \code{Tend} after \code{nGen} steps. 	
	T0			##<< the initial temperature (before the first step at iGen=0)
	, Tend=1	##<< the temperature at the last step
	, nGen		##<< the number of genrations	
	, iGen=1:nGen ##<< the steps for which to calculate the Temperature	
){
	# calcDEMCTemp
	##seealso<< 
	## \code{\link{calcTemperatedLogDen.default}}
	## \code{\link{twDEMC}}
	if( nGen < 1) return( numeric(0) )
	b = T0
	a = log(Tend/T0)/nGen
	### vector of Temperatures corresponding to steps iGen
	b*exp( a*iGen )
}
attr(calcDEMCTemp,"ex") <- function(){
	plot( 1:100, calcDEMCTemp(T0=100,Tend=5,nGen=100) )	
}

calcComponentTemp <- function(
	### calculating the temperature of result components of logDensity
	temp	##<< numeric scalar >= 1: global temperature
	,TFix	##<< named vector: temperature for the components that does not change but is held fixed
	,TProp	##<< named numeric vector [0,1]: temperature proportions of result components determines names and lenght of result
	,useMultiT=TRUE	##<< if set to FALSE only only global temperatue and TFix are applied
	,posTFix=match(names(TFix),names(TProp)) ##<< position index of TFix in vector of temperatures, specifiy for performance
){
	Ti <- if( useMultiT){
			structure( pmax(1, temp*TProp), names=names(TProp))
		}else{
			structure(rep(temp,length(TProp)),names=names(TProp))
		}
	if( 0 < length(posTFix)){ Ti[posTFix] <- TFix	}
	Ti
}


calcDEMCTempProp <- function(
	### Calculate Temperature of components 
	temp	##<< the maximum temperature
	,diffLogDen		##<< expected difference in LogDensitys proposed-accepted per datastream 
	,rFracMin=1/4	##<< fraction of max DiffDensity  below which temperatue is scaled down to yield  larger importance
){
	#rr <- diffLogDen/max(min(diffLogDen),1e-8)	# diffLogDen per largest misfit (lowest neg diff-LogDensity)
	#rr <- diffLogDen/min(diffLogDen)	# diffLogDen per largest misfit (lowest neg diff-LogDensity)
	rr <- max(diffLogDen/min(diffLogDen), 1e-8)	# ratio of diffLogLik of component per largest abs(diffLogLik), heed that all values are negative, hence min, ratio not below 1e-8
	#Ti <- pmin(temp,1+(temp-1)/rFracMin*rr)
	Ti <- pmin(temp,1+(temp-1)*(rr/rFracMin))
	Ti[rr<=0] <- 1	#give NA or negative values for Ti		
	Ti
	### vector of temperatures corresponding to diffLogDen with maximum corresponding to temp
}


