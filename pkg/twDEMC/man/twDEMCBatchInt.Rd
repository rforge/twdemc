\name{twDEMCBatchInt}
\alias{twDEMCBatchInt}
\title{twDEMCBatchInt}
\description{Calls \code{\link{twDEMCInt}} successively with \code{nBatch} generations.}
\usage{
twDEMCBatchInt(Zinit, nGen = 10, nBatch = 512, ..., restartFilename = NULL, fCheckConvergence = function(res, addArgs) {
    FALSE
}, fCheckConvergenceArgs = list(), T0 = 1, nGenBurnin = 0, doResetOutlierN = 256, nPops = 1, probUpDirBurnin = 0.75, controlTwDEMC = list(), doRepeatLowAcceptanceChains = TRUE, maxNGenBurnin = 50000, fCalcTStreamDiffLogLik = calcDEMCTempDiffLogLik3, argsFCalcTStreamDiffLogLik = list(), fCalcTGlobal = calcDEMCTempGlobal2, argsFCalcTGlobal = list())
}
\arguments{
  \item{Zinit}{the initial population}
  \item{nGen}{number of generations in total}
  \item{nBatch}{number of generations between saving results}
  \item{\dots}{further arguments passed to twDEMC}
  \item{restartFilename}{name of the file to save intermediate results resRestart.twDEMC, NULL for not saving}
  \item{fCheckConvergence}{checking convergence of a DEMC and interrupting }
  \item{fCheckConvergenceArgs}{additional arguments to the DEMC convergence checking function}
  \item{T0}{initial temperature of burnin , defaults to 1 or if Zinit is twDEMC to the temperature of the last row}
  \item{nGenBurnin}{initial number of generations of burnin (Temperature decrease to 0, and usage of probUpDirBurnin)
either numeric vector of length nPops or scalar}
  \item{doResetOutlierN}{if > 0, outlier chains are reset to best chain}
  \item{nPops}{number of independent populations}
  \item{probUpDirBurnin}{probUbDir during burnin (see twDEMC argument propuUpDir)}
  \item{controlTwDEMC}{controls to twDEMC, some items are overwritten}
  \item{doRepeatLowAcceptanceChains}{
}
  \item{maxNGenBurnin}{maximum burnin beyond which can not be extendend on too low acceptance rate}
  \item{fCalcTStreamDiffLogLik}{function calculate optimal Temperature and acceptance rates based on Lp-La}
  \item{argsFCalcTStreamDiffLogLik}{further arguments to fCalcTStreamDiffLogLik}
  \item{fCalcTGlobal}{function calculating global target temperature}
  \item{argsFCalcTGlobal}{further arguments to fCalcTStreamDiffLogLik}
}
\details{Usually invoked by \code{\link{twDEMCBatch}}.

If Zinit is of class twDEMC, initial temperature is set to the temperature of the last row
and the number of generations already in Zinit are skipped.

\describe{\item{Temperature estimate from proposal distribution}{
The distribution of differences between Likelihood of proposals Lp and of accepted state La
can be used to estimate an optimal temperature per data stream, so that each
datastream contributes to rejections in about the same magnitude and the overall
acceptance rate is about a specified value.
The proportions of the so calculated datastream specific temperature are multiplied 
with the global temperature on Metropolis decisions.
deprecated[Further, if the temperatures of the datasteams are all below the goal of the global
temperature Tend, Tend is also lowered]
}}

\describe{\item{Saving and restarting}{ 
If \code{restartFilename} has been specified then result are stored as variable \code{resRestart.twDEMC} to the given file.
Runs can be continued without the need of respecifying all the parameters (see example). 
This is helpful with expensive models and long cluster runs for cases where the program has to be aborted.
}}

\describe{\item{cooling and acceptance rate}{ 
If acceptance rate of some population drops below rate=minAccepRateTempDecrease then cooling is too fast.
In this moderate case do not repeat the rund but keep the current temperature for the next period for this population.
and extend the burnin phase by the length of this period.
}}}
\value{List of class twDEMC (see \code{\link{twDEMCInt}}) with additional entry nGenBurnin}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twDEMCInt}}
\code{\link{twDEMCBatch}}}
\examples{
data(twLinreg1); attach( twLinreg1 ) 

# run with saving a restart file each 50 generations
.nPops=2
argsFLogLik <- list(
	fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
	obs=obs,				### vector of data to compare with
	invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
	thetaPrior= thetaTrue,	### the prior estimate of the parameters
	invCovarTheta = invCovarTheta,	### the inverse of the Covariance of the prior parameter estimates
	xval=xval
)
Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=4*.nPops, nPops=.nPops)
restartFilename="exampleTwDEMCBatch_saveAndRestart.RData"
unlink(restartFilename)
res <-  twDEMCBatch( Zinit=Zinit, nGen=60,	nPops=.nPops
	, fLogLik=logLikGaussian, argsFLogLik=argsFLogLik
	, nBatch=50, restartFilename=restartFilename	# save each 50 generations
)

# load the restart file and continue
load( file=restartFilename )	# variable resRestart.twDEMC
calcNGen(resRestart.twDEMC)		# 48, last thinnging interval (each thin=4) before 50
res2 <- twDEMCBatch(resRestart.twDEMC)	# continue without needing to respecify parameters
calcNGen(res2)					# 60 as nGen
res3 <- twDEMCBatch(res2, nGen=100)	    # continue even further without needing to respecify parameters
calcNGen(res3)					# 100 as nGen

detach()
}
