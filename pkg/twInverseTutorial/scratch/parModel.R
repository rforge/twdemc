kinrespModel1 <- function(
	### Respiration for given microbial parameters at given time.
	param   ##<< parameters: named numeric vector(x0,r0,mumax)
	,time	##<< time (numeric vector)
	,lambda=0.9	##<< Ratio of growth associated (coupled) specific respiration to total specific respiration. Usually 0.9.
	,YCO2=1.5	##<< Ratio of assimilated carbon per respired carbon.  Usually 1.5.
){
	##details<< components of param: \describe{
	## \item{x0}{ initial microbial biomass (numeric scalar)
	## \item{r0}{ initial microbial activity (numeric scalar)
	## \item{mumax}{ maximum growth rate (numeric scalar)
	##}
	# resp <- x0*(1-r0)*(1/lambda-1)*mumax/YCO2 + x0*r0*1/lambda*mumax/YCO2 * exp(mumax*time)
	resp <- param["x0"]*(1-param["r0"])*(1/lambda-1)*param["mumax"]/YCO2 + param["x0"]*param["r0"]*1/lambda*param["mumax"]/YCO2 * exp( param["mumax"] * time)
	### respiration at given time points (numeric vector)
}

.tmp.f <- function(){
	ds2 <- subset(ds, time<20)
	fitLS6 <- gnls( resp ~ kinrespModel1(x0,r0,mumax,time), ds2
		,params=x0+r0+mumax~1
		,start=c(x0=140, r0=0.1, mumax=0.24)
		,weights=varPower(fixed=0.0)
	)
	coef(fitLS6)
	coef(fitLS)			# compare
	plot((resid(fitLS6,type="pearson"))~ds2$time)	# seems ok for weighted residuals
	
	plot( resp~time, data=ds2)
	lines(fitted(fitLS6)~ds2$time)
}

