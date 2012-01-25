# slowly decreasing Temperature (cost reduction factor)

twDEMCSA <- function(
	### simulated annealing DEMC 
	thetaPrior			##<< vector of parameters, point estimate
	,covarTheta			##<< the a prior covariance of parameters 
	,nChainPop=4		##<< number of chains within population
	,nPop=2				##<< number of populations
	,doIncludePrior=TRUE	##<< should the prior be part of initial population

	,dInfos				##<< argument to \code{\link{twDEMCBlockInt}}
	,qTempInit=c(0.4,0.2)	##<< quantile of logDensities used to calculate initial Temperature
	,nGen0=1024			##<< number of generations in the initial batch
	,... 				##<< further argument to \code{\link{twDEMCBlockInt}}
){
	#mtrace(initZtwDEMCNormal)
	Zinit0 <- initZtwDEMCNormal( thetaPrior, covarTheta, nChainPop=nChainPop, nPop=nPop, doIncludePrior=doIncludePrior)
	ss <- stackChains(Zinit0)
	
	logDenL <- lapply( dInfos, function(dInfo){
		resLogDen <- twCalcLogDenPar( dInfo$fLogDen, ss, argsFLogDen=dInfo$argsFLogDen)$logDenComp
	})
	# replace missing cases
	logDenDS <- abind( logDenL )
	boFinite <- apply(is.finite(logDenDS), 1, all )
	m0FiniteFac <- sum(boFinite) / nrow(logDenDS)
	if( m0FiniteFac < 0.9 ){
		## XXTODO extend Zinit0
	}
	Zinit <- Zinit0
	
	#------ intial temperatures: 0.4 percentile of rLogDen
	temp0 <- apply( logDenDS, 2, function(rLogDen){
		pmax(1, quantile( (max(rLogDen) - rLogDen),c(qTempInit) ))
	})

	res2 <- res <- twDEMCBlock( Zinit
		, nGen=nGen0
		#,nGen=64, debugSequential=TRUE
		, dInfos=dInfos
		, T0=temp[1]
		, TEnd=temps[3]
		, nPop=.nPop
	)

	
	
}
attr(twDEMCSA,"ex") <- function(){
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	
	dInfos=list(
		dSparce=list(fLogDen=denSparcePrior, argsFLogDen=list(thresholdCovar=0, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=0,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		a=list(dInfoPos="dSparce", compPos="a")
		,b=list(dInfoPos="dRich", compPos="b")
	)
	
	do.call( dInfos$dSparce$fLogDen, c(list(theta=twLinreg1$theta0),dInfos$dSparce$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twLinreg1$theta0),dInfos$dRich$argsFLogDen))
}

