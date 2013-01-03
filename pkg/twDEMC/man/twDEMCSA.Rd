\name{twDEMCSA}
\alias{twDEMCSA}
\title{twDEMCSA}
\description{simulated annealing DEMC }
\usage{twDEMCSA(thetaPrior, covarTheta, nGen = 512, nObs, ..., dInfos, 
    m0 = calcM0twDEMC(length(thetaPrior), nChainPop), controlTwDEMC = list(), 
    debugSequential = FALSE, restartFilename = NULL, nChainPop = 4, 
    nPop = 2, doIncludePrior = FALSE, ctrlBatch = list(nGen0 = m0 * 
        controlTwDEMC$thin * 3, useSubspaceAdaptation = FALSE), 
    ctrlT = list(TFix = numeric(0), TMax = numeric(0), qTempInit = 0.4), 
    ctrlConvergence = list(), ctrlSubspaces = list())}
\arguments{
  \item{thetaPrior}{vector of parameters, point estimate}
  \item{covarTheta}{the a prior covariance of parameters}
  \item{nGen}{number of generations in the initial batch, default 512}
  \item{nObs}{integer vector (nResComp) specifying the number of observations for each result component}
  \item{\dots}{further argument to \code{\link{twDEMCBlockInt}}}
  \item{dInfos}{argument to \code{\link{twDEMCBlockInt}}}
  \item{m0}{minimum number of samples in step for extending runs}
  \item{controlTwDEMC}{list argument to \code{\link{twDEMCBlockInt}} containing entry thin}
  \item{debugSequential}{set to TRUE to avoid parallel execution, good for debugging}
  \item{restartFilename}{filename to write intermediate results to}
  \item{nChainPop}{number of chains within population}
  \item{nPop}{number of populations}
  \item{doIncludePrior}{should the prior be part of initial population
<< Recommendation to set to false, because if the TRUE parameter is in initial set, the Temperature is set to 1}
  \item{ctrlBatch}{list of arguments controlling batch executions, see \code{\link{twDEMCSACont}} and \code{\link{divideTwDEMCSACont}}\describe{
\item{nGen0}{number of generations for the intial run}
\item{useSubspaceAdaptation}{if TRUE then overall space is devided and each subspace is explored with locally adapted DEMC, see	\code{\link{divideTwDEMCSACont}}}
}}
  \item{ctrlT}{list of arguments controlling Temperature decrease, see \code{\link{twDEMCSACont}}  and \code{\link{divideTwDEMCSACont}}}
  \item{ctrlConvergence}{list or arguments controlling check for convergence, see \code{\link{twDEMCSACont}}  and \code{\link{divideTwDEMCSACont}}}
  \item{ctrlSubspaces}{list of arguments controlling splitting and merging of subspaces, see \code{\link{divideTwDEMCSACont}}}
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
res <- res0 <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs
	, nGen=3*256
	, ctrlT=list( TFix=c(parmsSparce=1) )
	, ctrlBatch=list( nGenBatch=256 )
	, debugSequential=TRUE
	#, restartFilename=file.path("tmp","example_twDEMCSA.RData")
)
res <- twDEMCSA( res0, nGen=4*256 )	# extend the former run

(TCurr <- getCurrentTemp(res))
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
