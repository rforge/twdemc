data(delta14Catm)
c14Constants <- list(
	#iR14CUnit = 1e-12
	#,iR14CStandard = 1.176	# http://www.nosams.whoi.edu/clients/data.html,  (Karlen, et. al., 1968),
	iR14CUnit = 1.176*1e-12		##<< atomic ratio of standard, i.e. output in fraction modern
	,iR14CStandard = 1			# http://www.nosams.whoi.edu/clients/data.html,  (Karlen, et. al., 1968), 
	,yr14CStandard=1950
	,lambda=1/8267	#radioactive decay constant
	,delta14Catm=delta14Catm
)

calcLagged14CSeries <- function(
	### Calculate the 14C flux by applying the atmospheric 14C concentration with a time lag.
	inputYr			##<< the year of the carbon flux 
	,inputValue		##<< the value of the carbon flux
	,delta14Catm	##<< dataframe columsn 1: yr and 2: delta14C signal of the atmosphere
	,lag=0			##<< integer scalar: time lage in yr
	,iR14CStandard=c14Constants$iR14CStandard ##<< isotopic ratio of the oxalic acid standard (heed that units of isotopes differ to scale to same magnitude) 
	,delta13C=-25	##<< isotopic ratio of 13C in delta units for correction for fractionation
){
	#shift atmospheric curve, repeat first value for initial gap
	lagged14Catm <- delta14Catm
	if( lag > 0 ){
		lagged14Catm[,"delta14C"] <- c( rep(lagged14Catm[1,"delta14C"],lag), lagged14Catm[1:(nrow(delta14Catm)-lag),"delta14C"]) 
	}
	yrs <- inputYr
	delta13Ccorrected14C(inputValue * delta2iR14C( lagged14Catm[lagged14Catm[,"yr"] %in% yrs,"delta14C"], iR14CStandard=iR14CStandard ))
	### atomic ratio of input as multiple of the standard (dimensionless) corrected for fractionation.
} 

delta13Ccorrected14C <- function(
	### correct 14C atomic ratio for fractionation by delta13C value 
	iR				##<< 14C isotopic ratio, or fraction modern
	,delta13C=-25	##<< isotopic ratio of 13C in delta units for correction for fractionation
){
	if( delta13C == -25) iR else iR*((1-25/1000)/(1+delta13C/1000))^2 
}

delta2iR14C <- function(
	### Convert delta unit to atomic ratio.
	delta	##<< numeric vector of delta values
	,iR14CStandard=c14Constants$iR14CStandard	##<< numeric scalar atomic ratio of the standard
){
	(delta/1000+1)*iR14CStandard	#isotopic ratio not corrected for fractionation
	
}

ir14C2delta <- function(
	### Convert atomic ratio to delta units.
	iR	##<< numeric vector of atomic ratios
	,iR14CStandard=c14Constants$iR14CStandard
){
	(iR/iR14CStandard-1)*1000
}

delta2FM <- function(
	### Convert delta unit to fraction modern.
	delta	##<< numeric vector of delta values
){
	(delta/1000+1)	
}

FMC2delta <- function(
	### Convert fraction modern to delta units.
	fm	##<< numeric vector of atomic ratios
){
	(fm-1)*1000
}



decayIR14C <- function(
	### Calculate the atomic ratio for given years assuming only radioactive decay.
	yr		##<< numeric vector for years
	, iR0=c14Constants$iR14CStandard	##<< isotopic ratio at yr0
	, yr0=c14Constants$yr14CStandard	##<< year when iR0 is given
	, lambda=c14Constants$lambda		##<< decay constant
){
	iR0*exp(lambda*(yr0-yr))
}
attr(decayIR14C,"ex") <- function(){
	yr=seq(-10000,2000,length.out=30)
	plot( decayIR14C(yr)~yr )
}