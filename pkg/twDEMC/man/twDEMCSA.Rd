\name{twDEMCSA}
\alias{twDEMCSA}
\title{twDEMCSA}
\description{simulated annealing DEMC }
\usage{twDEMCSA(thetaPrior, covarTheta, nChainPop = 4, nPop = 2, doIncludePrior = FALSE, 
    dInfos, qTempInit = 0.4, nObs, nGen = 512, TFix = numeric(0), 
    TMax = numeric(0), nBatch = 4, maxRelTChange = 0.025, maxLogDenDrift = 0.3, 
    restartFilename = NULL, TDecProp = 0.8, critSpecVarRatio = 20, 
    ...)}
\arguments{
  \item{thetaPrior}{vector of parameters, point estimate}
  \item{covarTheta}{the a prior covariance of parameters}
  \item{nChainPop}{number of chains within population}
  \item{nPop}{number of populations}
  \item{doIncludePrior}{should the prior be part of initial population
<< Recommendation to set to false, because if the TRUE parameter is in initial set, the Temperature is set to 1}
  \item{dInfos}{argument to \code{\link{twDEMCBlockInt}}}
  \item{qTempInit}{quantile of logDensities used to calculate initial beginning and end temperature, with default 0.4: 40\% of the space is accepted}
  \item{nObs}{integer vector (nResComp) specifying the number of observations for each result component}
  \item{nGen}{number of generations in the initial batch}
  \item{TFix}{numeric vector (nResComp) specifying a finite value for components with fixed Temperatue, and a non-finite for others
<< Alternatively a named numeric vector specifying the fixed temperature only for certain components. Note this breaks if result components (possibly of different logDensity functions) have same names}
  \item{TMax}{numeric vector (nResComp) specifying a maximum temperature for result components.
<< Alternatively a named numeric vector specifying the maximum temperature only for certain components. Note this breaks if result components (possibly of different logDensity functions) have same names}
  \item{nBatch}{number of batches with recalculated Temperature}
  \item{maxRelTChange}{if Temperature of the components changes less than specified value, we do not need further batches because of Temperature}
  \item{maxLogDenDrift}{if difference between mean logDensity of first and fourth quartile of the sample is less than this value, we do not need further batches because of drift in logDensity}
  \item{restartFilename}{filename to write intermediate results to}
  \item{TDecProp}{proportion of Temperature decrease: below one to diminish risk of decreasing Temperature too fast (below what is supported by other data streams)}
  \item{critSpecVarRatio}{if proprotion of spectral Density to Variation is higher than this value, signal problems and resort to subspaces}
  \item{\dots}{further argument to \code{\link{twDEMCSA}}}
}
\details{\describe{\item{initial temperature}{
Initial parameters are ranked according to their maximum log-Density across components.
The parameters at rank positions defined by qTempInit are selected.
The minimum log-Density or each datastream is used obtain initial temperature as: T =  -2/nObs logDen
}}}


\author{Thomas Wutzler}




\examples{
data(twTwoDenEx1)

thetaPrior <- twTwoDenEx1$thetaTrue
covarTheta <- diag((thetaPrior*0.3)^2)
invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation

thresholdCovar = 0.3	# the true value
thresholdCovar = 0		# the effective model that glosses over this threshold

dInfos=list(
	dSparce=list(fLogDen=denSparcePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
	,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
)
blocks = list(
	a=list(dInfoPos="dSparce", compPos="a")
	,b=list(dInfoPos="dRich", compPos="b")
)

do.call( dInfos$dSparce$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparce$argsFLogDen))
do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))

#str(twTwoDenEx1)
nObs <- c( parmsSparce=1, obsSparce=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )

#trace(twDEMCSACont, recover )
#trace(twDEMCSA, recover )
resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs
	, TFix=c(parmsSparce=1)
	, nGen=256
	#, nBatch=2
	, nBatch=8
	, debugSequential=TRUE
	#, restartFilename=file.path("tmp","example_twDEMCSA.RData")
)

(TCurr <- getCurrentTemp(resPops))
mc0 <- concatPops(res)
plot( as.mcmc.list(mc0) , smooth=FALSE )
matplot( mc0$temp, type="l" )
logDenT <- calcTemperatedLogDen(stackChains(mc0$resLogDen), TCurr)
iBest <- getBestModelIndex( logDenT, res$dInfos )
maxLogDenT <- logDenT[iBest, ]
ss <- stackChains(mc0$parms)
(thetaBest <- ss[iBest, ])
(.qq <- apply(ss,2,quantile, probs=c(0.025,0.5,0.975) ))
plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(rowSums(logDenT),c(10,100))] )
plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"obsSparce"],c(10,100))] )
plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"logDen1"],c(10,100))] )
plot( ss[,"a"], ss[,"b"], col=rgb(
		twRescale(logDenT[,"obsSparce"]),0, twRescale(logDenT[,"logDen1"]) ))
apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )

# density of parameters
plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")

# predictive posterior (best model only)
pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 

}
