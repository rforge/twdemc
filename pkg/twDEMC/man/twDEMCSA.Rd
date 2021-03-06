\name{twDEMCSA}
\alias{twDEMCSA}
\title{twDEMCSA}
\description{simulated annealing DEMC }
\usage{twDEMCSA(thetaPrior, covarTheta, nGen = 512, nObs, ..., dInfos, 
    m0 = calcM0twDEMC(length(thetaPrior), nChainPop), controlTwDEMC = list(), 
    debugSequential = FALSE, restartFilename = NULL, remoteDumpfileBasename = NULL, 
    nChainPop = 4, nPop = 2, doIncludePrior = FALSE, ctrlBatch = list(nGen0 = m0 * 
        controlTwDEMC$thin * 3), ctrlT = list(qTempInit = 0.4, 
        TBaseInit = NULL, isVerbose = FALSE, TEndFixed = NULL, 
        qBest = 0.1), ctrlConvergence = list())}\arguments{
  \item{thetaPrior}{vector of parameters, point estimate
<< , alternatively array with initial states, as returned by \code{\link{initZtwDEMCNormal}}}
  \item{covarTheta}{the a prior covariance of parameters, see \code{link{initZtwDEMCNormal}}}
  \item{nGen}{number of generations}
  \item{nObs}{integer vector (names resComp) specifying the number of observations for each result component
for each resCompName there must be an entry in nObs}
  \item{\dots}{further argument to \code{\link{twDEMCBlockInt}}}
  \item{dInfos}{argument to \code{\link{twDEMCBlockInt}}}
  \item{m0}{minimum number of samples in step for extending runs

minimum number of samples in step for extending runs}
  \item{controlTwDEMC}{list argument to \code{\link{twDEMCBlockInt}} containing entry thin}
  \item{debugSequential}{set to TRUE to avoid parallel execution, good for debugging}
  \item{restartFilename}{filename to write intermediate results to}
  \item{remoteDumpfileBasename}{fileBasename to write dumps to on error}
  \item{nChainPop}{number of chains within population}
  \item{nPop}{number of populations}
  \item{doIncludePrior}{should the prior be part of initial population
<< Recommendation to set to false, because if the TRUE parameter is in initial set, the Temperature is set to 1}
  \item{ctrlBatch}{list of arguments controlling batch executions, see \code{\link{twDEMCSACont}}\describe{
\item{nGen0}{number of generations for the initial batch}
}}
  \item{ctrlT}{list of arguments controlling Temperature decrease, see \code{\link{twDEMCSACont}}\describe{
\item{qTempInit}{quantile of logDensities used to calculate initial beginning and end temperature, with default 0.4: 40\% of the space is accepted}
\item{TBaseInit}{numeric scalar: initial base temperature. If given this is used for calculating initial temperature}
\item{isVerbose}{boolean scalar: set to TRUE to report stream base temperatures during batches}
\item{TEndFixed}{set to a scalar end temperature, e.g. in order to decrease temperatue to a given temperature}
\item{qBest}{the proportion of best samples to base calculation of target temperature on}
}}
  \item{ctrlConvergence}{list or arguments controlling check for convergence, see \code{\link{twDEMCSACont}}}
}
\details{\describe{\item{Initial temperature}{
Initial parameters are ranked according to their maximum log-Density across components.
The parameters and logDensity results at rank position defined by argument \code{ctrlT$qTempInit} is selected.
The stream temperatures are inferred by deviding logDensity components by the number of observations.
From these, a common base temperatue is calculated by (see \code{\link{calcBaseTemp}}) and rescaled to stream temperatures by \code{\link{calcStreamTemp}}.
}}}
\value{An object of class \code{twDEMCPops} as described in  \code{\link{twDEMCBlockInt}}.
See \code{\link{subset.twDEMCPops}} for processing further handling of this class.}

\author{Thomas Wutzler}




\examples{
    
#--------------- single density ----------------------------
# we will use logDenGaussian as logDensity function that compares a simple linear model with observations
data(twLinreg1)

# collect all the arguments to the logDensity in a list (except the first argument of changing parameters)    
argsFLogDen <- with( twLinreg1, list(
        fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
        obs=obs,			    ### vector of data to compare with
        invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
        thetaPrior = thetaTrue,	### the prior estimate of the parameters
        invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
        xval=xval
))
do.call( logDenGaussian, c(list(theta=twLinreg1$theta0),argsFLogDen))
do.call( logDenGaussian, c(list(theta=twLinreg1$thetaTrue),argsFLogDen))    # slightly largere misfit than nObs/2=15, underestimated sdObs

.nGen=200
.nPop=2
mcPops <-  twDEMCSA( 
        theta=twLinreg1$theta0, covarTheta=diag(twLinreg1$sdTheta^2)       # to generate initial population
        , nGen=.nGen
        , dInfos=list(den1=list(fLogDen=logDenGaussian, argsFLogDen=argsFLogDen))
        , nPop=.nPop                                        # number of independent populations
        , controlTwDEMC=list(thin=4)                        # see twDEMCBlockInt for all the tuning options
        , ctrlConvergence=list(maxRelTChangeCrit=0.1)       # ok if T changes less than 10\% 
        , ctrlT=list(TFix=c(parms=1))                       # do not use increased temperature for priors
        , nObs=c(obs=length(argsFLogDen$obs))               # number of observations used in temperature calculation
)
#mcp <- twDEMCSA( mcp, nGen=2000) 
mcPops <- twDEMCSA( mcPops, nGen=400)     # continue run 

rescoda <- as.mcmc.list(mcPops)
plot(rescoda, smooth=FALSE)
mcChains1 <- concatPops(mcPops)                   # array representation instead of list of pops, last dim is the chain
mcChains2 <- concatPops(stackChainsPop(mcPops))   # combining dependent chains within one population
mcChains3 <- concatPops(subsetTail(mcPops,0.5))   # take only the last part of the chains
c(getNGen(mcChains1), getNGen(mcChains3))
plot(as.mcmc.list(mcChains2), smooth=FALSE)



#--------------- multiple densities -------------------------
data(twTwoDenEx1)

thetaPrior <- twTwoDenEx1$thetaTrue
covarTheta <- diag((thetaPrior*0.3)^2)
invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation

thresholdCovar = 0.3	# the true value used to generate the observations
thresholdCovar = 0		# the effective model that glosses over this threshold

#str(twTwoDenEx1)
nObs <- c( obsSparse=length(twTwoDenEx1$obs$y1), obsRich=length(twTwoDenEx1$obs$y2) )

dInfos=list(
	dSparse=list(fLogDen=denSparsePrior
        , argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta)
        #, maxLogDen=-1/2*nObs[c("parmsSparse","obsSparse")] # control overfitting
        )
	,dRich=list(fLogDen=denRichPrior
        , argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta)
        , maxLogDen=c(parmsRich=0,-1/2*nObs["obsRich"])     # control overfitting for rich datastream
        )
)
blocks = list(
	a=list(dInfoPos="dSparse", compPos="a")
	,b=list(dInfoPos="dRich", compPos="b")
)

names(do.call( dInfos$dSparse$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparse$argsFLogDen)))
names(do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen)))


#trace(twDEMCSACont, recover )
#trace(twDEMCSA, recover )
res <- res0 <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs
	, nGen=3*256
	, ctrlT=list( TFix=c(parmsSparse=1,parmsRich=1) )   # no increased Temperature for priors
	, ctrlBatch=list( nGenBatch=256 )
	, debugSequential=TRUE
    , controlTwDEMC = list(
           DRgamma=0.1                          # use Delayed rejection
           #,controlOverfittingMinNObs = 20      # use overfitting control (for obsRich), recommended on using single density 
    )
	#, restartFilename=file.path("tmp","example_twDEMCSA.RData")
)
res <- twDEMCSA( res0, nGen=2*256 )	# extend the former run

(TCurr <- getCurrentTemp(res))
# Note that T does decrease to 1
# This accounts for structural model mismatch in addition to observation uncertinaty

mc0 <- concatPops(res)
mcE <- concatPops(subsetTail(res,0.2))      # only the last 20\%
plot( as.mcmc.list(mc0) , smooth=FALSE )
matplot( mc0$temp, type="l" )
logDenT <- calcTemperatedLogDen(stackChains(mcE$resLogDen), TCurr)
iBest <- getBestModelIndex( logDenT, res$dInfos )
maxLogDenT <- logDenT[iBest, ]
ss <- stackChains(mcE$parms)
(thetaBest <- ss[iBest, ])
twTwoDenEx1$thetaTrue
(.qq <- apply(ss,2,quantile, probs=c(0.025,0.5,0.975) ))    # model error really in paramter b
plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twMisc::twRescale(rowSums(logDenT),c(10,100))] )
plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twMisc::twRescale(logDenT[,"obsSparse"],c(10,100))] )
plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twMisc::twRescale(logDenT[,"obsRich"],c(10,100))] )
plot( ss[,"a"], ss[,"b"], col=rgb(
		twMisc::twRescale(logDenT[,"obsSparse"]),0, twMisc::twRescale(logDenT[,"obsRich"]) ))
apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )

# density of parameters
plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")

# predictive posterior (best model only)
pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparse=xSparse, xRich=xRich) )
plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1)    # note that deviation is now in y2 - consistent with the introduced bias
}
