initState.howland.ICBM1SteadyState  <- function(
	### Initial states from parms for ICBM1 in steady state
	padj		##<< must contain cY, kO, Ctot, yr0
	,modMeta=modMetaICBM1()
	,delta14Catm=c14Constants$delta14Catm
){
	iRNew <- delta2iR14C(delta14Catm$delta14C[delta14Catm$yr==padj$yr0])
	tvrOld <- 1/padj$kO	#1000 yr old carbon
	iROld <- decayIR14C( yr=padj$yr0, iR0=delta2iR14C(delta14Catm$delta14C[1]), yr0=1950-tvrOld )	# near 1 (standard of old wood)
	#mtrace(initStateSoilMod)
	x0 <- initStateICBM1( xc12=padj$Ctot*c(padj$cY,(1-padj$cY)),iR=matrix(c(iRNew,iROld),ncol=1,dimnames=list(NULL,"c14")) )
}

meanInput <- function(
	### provide and unperturbed data series, neglect the standard deviation
	input	##<< list of datastream matrices with first two columns time and value 
	,padj	##<< parameters
){
	lapply(input, "[", ,1:2, drop=FALSE)
}

of.howlandSteady <- function(
	### Objective function for comparing against Howland data, assuming constant input and steady state C-stocks, which determines k-Values.
	normpopt	##<< numeric vector: the point in normal parameter space
	, logLikAccept=numeric(0)	##<< numeric named vector the logLikelihood of previously accepted parameters, see details
	, metropolisStepTemp=structure(rep(1.0,length(logLikAccept)),names=names(logLikAccept))
	### the temperature for current metropolis decision
	, ...			##<< further parameters to model$fSolve
	, model		
	### model to evaluate. A list with following components \describe{
	### \item{modMeta}{meta information about the model. see \code{\link{twCreateModMeta}} }
	### \item{fInitState}{function to initialize state based on parameters and amendment. see \code{\link{initStateSoilMod}} }
	### \item{fSolve}{function to ODE across time. see \code{\link{solveICBM1}} }
	### }
	, poptDistr	
	### information on parameter distributions, list with entries \describe{
	### \item{trans}{character vector: type of distribtution (norm,lognorm,logitnorm)}
	### \item{mu}{numeric vector: distribution parameter mu, i.e. expected values at normal scale}
	### \item{invsigma}{numeric matrix: inverse of the distribution parameter sigma on multivariate (transformed) normal distribution}
	### }
	### See also \code{\link{twQuantiles2Coef}}.
	### Parameters order must correspond to normpopt. Use \code{\link{twConstrainPoptDistr}} to select a constrained subset from all parameters.
	, obs=Howland14C$obsNutrientSite	##<< see \code{\link{Howland14C}}
	#, times=c(1950, sort( unique(unlist(lapply(obs,"[",,1))) ))
	, times=1950:2007				##<< time points in yr to be modelled and compared.
	, input=Howland14C$litter		##<< carbon inputs to reservoir C see \code{\link{Howland14C}$litter}.
	, parms=HowlandParameterPriors$parms0		##<< default parameters (for the non-optimized ones)
	, fCalcBiasedObs=NULL			##<< function(obs,padj,...){obs} possibility to account for bias and to optimize bias parameters 
	, argsFCalcBiasedObs=list()		##<< further arguments to fCalcBiasedObs
	, fCalcBiasedInput=meanInput	##<< function(input,padj,...){obs} possibility to account for bias and to optimize bias parameters 
	, argsFCalcBiasedInput=list()	##<< further arguments to fCalcBiasedInput
	, popt.names=names(normpopt)	##<< names of the parameters. They are sometimes stripped by fitting algorithms.
	, fTransOrigPopt=transOrigPopt.default	##<< function that translates parameters from normal to original scale
	, doStopOnError=FALSE					##<< by default -Inf is returned on error, set to TRUE for debugging
	, includeStreams=c(names(obs),"parms") 	##<< character vector of subset of data streams to include
	, delta14Catm=c14Constants$delta14Catm	##<< 14C signature of the atmosphere: signature of litter inputs after a time lag.
	, fCalcIROLayer=calcIROLayer			##<< function to calculate iRofO-Layer
	#, facPrior=1							##<< factor to scale prior information (we have only weak data), detoriates the prior for transformed data so that DEMC does not work any more
){
	# ofb_FS.hamer
	##details<< objective function for fitting SoilMod_FS to respRate timeseries of Hamer incubation experiment of both control, amendm and c14obs
	## inclDataSeries can be set to a subset of 1:3 corresponding to control, amendment_large, and amendment_small respectively to save computing time
	
	#tmp <- matrix(1,3,3); tmp[4,]	#test dumping remote error
	if( is.null(popt.names) )
		stop("of.howlandSteady: normpopt provided without names attribute or argument popt.names")
	
	#datastreamNames <- c("parms", "control","controlSum","amendm","amendmSum","c14obs","c14obsSum","amendm_small","amendmSum_small","c14obs_small","c14obsSum_small") #names of the dataseries
	#each component for each observation with independent metropolis decisions
	#includeStreams=c("parms","control")
	nComp <- length(includeStreams)	#number of misfit components returned when sucessful
	misfit <- misfit0 <- structure( rep(NA_real_,nComp), names=includeStreams ) #initially set to 0, i.e. no constraint
	
	names(normpopt) <- popt.names
	popt = fTransOrigPopt(normpopt,poptDistr$trans)
	padj <- parms;	padj[popt.names] <- popt	# replace parameters to optimized in parms
	
	# single logLik components of parameters
	#parameter misfit
	diff.p <- normpopt - poptDistr$mu
	logLikParms <- structure( -0.5 * as.vector(diff.p^2 %*% poptDistr$invSigma), names=names(normpopt) )
	
	resSolve <- NULL
	misfitFail <- function(msg=NULL, errmsg=NULL){
		if( doStopOnError & !is.null(errmsg) ) stop(errmsg)
		if( !any(is.na(misfit)) ){
			#non-finite values cause rejection in .doDEMCStep
			#usually misfitFail is called before all components of misfit are changed away from NA.
			#if all components are initialized one is set to NA trigger rejection
			#in order to allow inspection of proposals avoid using the same component for this purpose
			misfit[ sample.int(length(misfit),1)] <- NA
		}
		#misfit[ extComp[sample.int(length(extComp),1)] ] <- Inf 	#will cause rejection
		rLogLik <- -1/2*misfit
		attr(rLogLik,"out") <- resSolve
		attr(rLogLik,"msg") <- msg
		attr(rLogLik,"errmsg") <- errmsg
		attr(rLogLik,"obs") <- obsadj
		attr(rLogLik,"logLikParms") <- logLikParms   
		return(rLogLik)
	}
	#mtrace(misfitFail)
	#check parameters outside range before calculating model
	if( any(!is.finite(popt)))	return( misfitFail(errmsg="Encountered nonfinite parameters.") )	#exponent function produced non-finite numbers
	
	if( "parms" %in% includeStreams) misfit["parms"] =  t(diff.p) %*% poptDistr$invSigma %*% diff.p
	
	##details<<
	## Supports a multi-step Metropolis descision. If \code{logLikAccept["parms"]} is provided, 
	## then a Metropolis descision is done based only on the parameters.
	## If it fails, then -Inf is returned. 
	## The possible costly evaluation of fModel is avoided.
	## \cr
	## Similarly, the small amendmend experiment and the control are not calculated, if 
	## observations are rejected be the large amendment experiment.
	## Provide a subset of \code{logLikAccept[c("amdend","c14obs","amdend_small","c14obs_small")]} to support early rejection.
	logrunif1 <- log(runif(1))  	#use the same random number of all metropolis decisions
	if( (!is.na(logLikXParms <- logLikAccept["parms"])) && (logLikXParms<0)){
		logLikPropParms = -1/2 * misfit["parms"]
		logr = (logLikPropParms - logLikXParms) / metropolisStepTemp["parms"]
		if ( is.numeric(logr) & (logr) <= logrunif1 )
			return( misfitFail("internal Metropolis rejection on parameters.") )
	}
	
	obsadj <- if( is.function(fCalcBiasedObs) ) do.call( fCalcBiasedObs, c(list(obs,padj),argsFCalcBiasedObs)) else obs	# biased observation data
	# assume steady state: litter inputs = respiration, recalculate root input
	input$root[1,"obs"] <- obsadj$respCum[1,"obs"] - input$leaf[1,"obs"]
	inputadj <- if( is.function(fCalcBiasedInput) ) do.call( fCalcBiasedInput, c(list(input,padj),argsFCalcBiasedInput)) else input	# biased input data
	
	#initial states for all three treatments
	padj$Ctot <- obs$somStock[1,2]
	sumInput <- sum(sapply(inputadj,"[",1,2))
	padj[c("kY","kO")] <- calcSteadyK_ICBM1(Ctot=padj$Ctot,cY=padj$cY,h=padj$h,iY=sumInput)
	
	padj$yr0 <- times[1]
	x0 <- model$fInitState(padj, modMeta=model$modMeta, delta14Catm=delta14Catm)
	
	checkModelErr <- function(resSolve){
		if( inherits(resSolve, "try-error")){
			misfit[1] <- -Inf	# NA indicates not accessed, if solving model fails this means Log-Likelihood of -Inf
			return(resSolve)
		} 
		if( any(is.na(resSolve[-1,])) ){	#omit first row because no cumulative resp, F14C NAN
			misfit[ 1 ] <- -Inf	# NA indicates not accessed, if solving model fails this means Log-Likelihood of -Inf
			return("solving model produced NaNs")
		} 
		if( (nrow(resSolve) != length(times))  ){
			misfit[ 1 ] <- -Inf	# NA indicates not accessed, if solving model fails this means Log-Likelihood of -Inf
			return(paste(
					"lsoda return different number of rows than time:"
					,", nrow(out0.tmp)=",nrow(resSolve),", length(times)=",length(times),sep=""))
		}
		character(0)
		### error msg if model output is invalid
		### character(0= otherwise
		### side effect: misfit[1] is set to -Inf
	}
	
	resSolve <- try( model$fSolve(x0=x0, times=times, parms=padj, input=inputadj, delta14Catm=delta14Catm, modMeta=model$modMeta, ...) )
	if( 0<length(errmsg<-checkModelErr(resSolve))) return( misfitFail(errmsg=errmsg))
	# estimate isotopic ratio of O-Layer
	resT <- resSolve	# calculate for all years
	iROLayer <- fCalcIROLayer(
		somStockMineral = obsadj$somStock[1,2]-obsadj$somOStock[1,2]
		,somStock = resT[,"cStock",drop=FALSE]
		,irS = resT[,"F14CT",drop=FALSE]
		,irM = resSolve[1,"F14C_O"]	#assume equals slow pool at beginning of the observation
	)
	if( "respCum" %in% includeStreams ){
		resT <- resSolve[ resSolve[,1] %in% obsadj$respCum[,"times"], ,drop=FALSE]	
		misfit["respCum"] <- sum( ((resT[,"R_c12"]/diff(range(times))-obsadj$respCum[,"obs"])/obsadj$respCum[,"sdObs"])^2 )
	}
	if( "respFM" %in% includeStreams ){
		resT <- resSolve[ resSolve[,1] %in% obsadj$respFM[,"times"], ,drop=FALSE]	
		misfit["respFM"] <- sum( ((resT[,"respF14CT"]-obsadj$respFM[,"obs"])/obsadj$respFM[,"sdObs"])^2 )
	}
	if( "somStock" %in% includeStreams ){
		resT <- resSolve[ resSolve[,1] %in% obsadj$somStock[,"times"], ,drop=FALSE]	
		misfit["somStock"] <- sum( ((resT[,"cStock",drop=FALSE]-obsadj$somStock[,"obs"])/obsadj$somStock[,"sdObs"])^2 )
	}
	if( "somOFM" %in% includeStreams ){
		#-- wrong assumption c14 old pool does not change: assumption mineral=old pool
		iROLayerObsYears <- iROLayer[ resSolve[,1] %in% obsadj$somOFM[,"times"], ,drop=FALSE]
		# add 10% uncertainty due to reconstruction of iRO
		misfit["somOFM"] <- sum( ((iROLayerObsYears-obsadj$somOFM[,"obs"])/(obsadj$somOFM[,"sdObs"]+padj$iROLayerCalcRelErr*obsadj$somOFM[,"obs"]))^2 )
		#-- does not work: assume old pool is mineral soil yong pools is organic layer
		#resT <- resSolve[ resSolve[,1] %in% obsadj$somOFM[,"times"], ,drop=FALSE]	
		#misfit["somOFM"] <- sum( ((resT[,"F14C_Y",drop=FALSE]-obsadj$somOFM[,"obs"])/obsadj$somOFM[,"sdObs"])^2 )
		#-- assumption: somOFM = somFM (same as fraction of total soil)
		#misfit["somOFM"] <- sum( ((resT[,"F14CT",drop=FALSE]-obsadj$somOFM[,"obs"])/obsadj$somOFM[,"sdObs"])^2 )
	}
	if( "somOStock" %in% includeStreams ){
		# not modelled
		misfit["somOStock"] <- 0
	}
	
	
	if( any(is.na(misfit))) stop("Encountered NA at the end of the of.howland.steadystate.")
	
	# weighted least squares, without the last observation (which is NA because of diff)
	#c(sum(t(diff.ctrl[-lc])^2),sum(t(diff.amendm[-llt])^2))
	rLogLik <- -1/2* misfit 	#twutz: adjusted 1/2 to comply to negative log-likelihood
	attr(rLogLik,"out") = resSolve
	attr(rLogLik,"parms") = padj	# includes calibrations of parameters (F0 to maximum etc)
	attr(rLogLik,"obs") = obsadj	# includes calibrations of parameters (F0 to maximum etc)
	attr(rLogLik,"input") = inputadj	# includes calibrations of parameters (F0 to maximum etc)
	attr(rLogLik,"iROLayer") = iROLayer 
	#t(diff.p) %*% poptDistr$invSigma %*% diff.p
	attr(rLogLik,"logLikParms") = logLikParms  
	#if( any(rLogLik == Inf) ) stop("InfError: rLogLik=",rLogLik)
	rLogLik
}


calcIROLayer <- function(
	### Calculate the isotopic ratio of the O-Layer.
	somStockMineral		##<< numeric scalar: the CStock of the mineral soil (assumed not to change)
	, somStock			##<< numeric vector: the CStock of the total soil
	, irS				##<< isotopic ratio of the total SOM carbon
	, irM				##<< isotopic ratio of the mineral SOM carbon (assumed not to change during this time)
){
	# mixing model irS = (1-cM)*irO + cM*irM)
	#    irO = (irS - cM*irM)/(1-cM)	
	# 	 cM: proportion of mineral soil to entire = (somStockMineral)/somStock
	#	 irM: isotopic ratio of mineral soil: assume initial of old pool
	#ifM <- x0
	cM <- somStockMineral/somStock
	irO <- (irS - cM*irM)/(1-cM)
	### numeric vector
}
attr(calcIROLayer,"ofCall") <- function(){
	resT <- resSolve	# calculate for all years
	iROLayer <- calcIROLayer(
		somStockMineral = obsadj$somStock[1,2]-obsadj$somOStock[1,2]
		,somStock = resT[,"cStock",drop=FALSE]
		,irS = resT[,"F14CT",drop=FALSE]
		,irM = resSolve[1,"F14C_O"]	#assume equals slow pool
	)
}

