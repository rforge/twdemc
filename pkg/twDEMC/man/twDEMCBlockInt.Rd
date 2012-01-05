\name{twDEMCBlockInt}
\alias{twDEMCBlockInt}
\title{twDEMCBlockInt}
\description{Differential Evolution Markov Chain with blocked parameter update}
\usage{
twDEMCBlockInt(pops = list(list(parms,
nGen = 12,
T0 = 1,
TEnd = 1,
TProp = NULL,
X = NULL,
logDenCompX = NULL,
spaceInd = 1,
upperParBounds = numeric(0),
lowerParBounds = numeric(0),
splits = numeric(0))),
dInfos = list(list(fLogDen = NULL,
xC = xC,
logDenCompC = logDenCompC,
parUpdateDenC = parUpdateDenC,
compPosDen = 1:nrow(pops[[1]]$parms),
argsFLogDen = list(),
intResComp = vector("integer",
0),
fLogDenScale = 1,
TFix = vector("numeric",
0))),
blocks = list(list(dInfoPos = 1,
compPos = dInfos[[1]]$compPosDen,
fUpdateBlock = updateBlockTwDEMC,
argsFUpdate = list(),
requiresUpdatedDen = TRUE)),
nGen = integer(0),
spacePop = 1:nPop,
controlTwDEMC = list(),
debugSequential = FALSE,
remoteDumpfileBasename = NULL,
fDiscrProp = NULL,
argsFDiscrProp = list(),
doRecordProposals = FALSE)
}
\arguments{
  \item{pops}{list of population infos for each population, each info is a list with components\describe{
\item{parms}{list of matrices (nState x nParm  x nChain) initial states for each population see details and \code{\link{initZtwDEMCNormal}}.}
\item{nGen}{scalar integer: number of generations, if given overwrites \code{pops[[i]]$nGen}}
\item{T0}{initial temperature}
\item{TEnd}{end temperature}
\item{TProp}{numeric vector (nResComp) temperature proportions of result components
can also be given as character vector with names of result components, however, be aware that htis fails if several logDen may return components of the same name}
\item{X}{numeric matrix (nParm x nChainPop) initial state}
\item{logDenCompX}{numeric matrix (nComp x nChain): logDen components of initial state X, see details}
\item{spaceInd}{the space replicate that this population belongs to}
\item{lowerParBounds}{similar to upperParBounds: sample > bound}
\item{splits}{named numeric vector of splitting points, used to remerge divided subspaces}
}}
  \item{dInfos}{named list of used density functions. Each entry is a list with components\describe{
\item{xC}{numeric vector: current accepted state}
\item{logDenCompC}{numeric vector: result components of fLogDen for current position}
\item{parUpdateDenC}{integer vector: logDensity that recently updated parameter at given index
to handle first steps in Multi-step Metropolis decision internally. 
See details.}
\item{compPosDen}{index or names of the parameter components that are used in this density function}
\item{argsFLogDen}{further arguments passed to fLogDen}
\item{TFix}{named numeric vector with Temperature for result components that have fixed temperature}
}}
  \item{blocks}{list of parameter blocks, each block is a list with entries\describe{
\item{dInfoPos}{name or position to \code{fLogDenInfo}. Several blocks may share the same density but update different parameters}
\item{compPos}{names or index of the parameter components to be updated}
\item{fUpdateBlock}{function to update the parameters.
\cr It must return a list with first three components xC, logDenCompC, and parUpdateDenC 
as described in \code{\link{updateBlockTwDEMC}}}
\item{argsFUpdate}{further arguments passed to fUpdate}
\item{requiresUpdatedDen}{if fUpdateBlock does not depend on current density result, then set this to FALSE and save some calculation time}
}}
  \item{nGen}{scalar integer: number of generations, if given overwrites \code{pops[[i]]$nGen}}
  \item{spacePop}{the space replicate that each population belongs to. By default assume only one population per space, overridden by entry in pops}
  \item{controlTwDEMC}{control parameters influencing the update and the convergens of the chains (see details)}
  \item{debugSequential}{if TRUE apply is used instead of sfApply, for easier debugging}
  \item{remoteDumpfileBasename}{the basename of a dumpfile that is created on error on remote process}
  \item{fDiscrProp}{function applied to proposal, e.g. to round proposals to to discrete possible values function(theta,...)}
  \item{argsFDiscrProp}{further arguments to fDiscrProp}
  \item{doRecordProposals}{if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y}
}
\details{This is the central method for applying a Differential Evolution Markov Chain given a function of 
an unnormalized probability density.
It is invoked usually by (\code{\link{twDEMCBlock.array}} or \code{\link{twDEMCBlock.twDEMCPops}})

\describe{ \item{recognized control parameters in \code{controlTwDEMC}: }{\describe{
\item{thin}{thinning interval}
\item{F}{related to multiplicative error (F2=F/sqrt(2*Npar), see eps.mult}
\item{pSnooker}{probability of a snooker update (others parallel updates)}
\item{pGamma1}{probability of jumping to state of another chain (different modes)}
\item{epsMult}{>0 gives d-dimensional variation around gamma. It adds scaled uncorrelated noise to the proposal. Its advantage over eps.add is that its effect scales with the differences of vectors in the population whereas eps.add does not. if the variance of a dimensions is close to 0, eps.mult gives smaller changes. \cr A uniformly distributed error, i.e. F2*runif(1+-epsMult*prop) multiplied to difference vector from parallel update}
\item{epsAdd}{>0 is needed to ensure that all positions in the space can be reached. For targets without gaps, it can set small or even to 0. \cr sd of normally distributed error added to proposal by parallel or snooker update.}
\item{pAcceptWindowWidth}{number of generations back over which the acceptance rate is calculated}
\item{probUpDir}{probability of direction between two states of increasing Density (increase during burin may accelerate convergence)}
\item{initialAcceptanceRate}{numeric matrix (nBlock x nChains) initially assumed acceptance rate. Used to calculate the number of generations backwards to sample from}
\item{DRgamma}{factor for reducing step length [0..1) in delayed rejection step, 0 means no DR step}
\item{minPCompAcceptTempDecr}{if acceptance rate drops below minPCompAcceptTempDecr+0.02 times this level, employ delayed rejection (DR)}
\item{pIndStep}{independent state about after on average about those number of 1.5 accepted steps}
\item{nPastGen}{factor for determining the number of recent past states to sample during burnin. It is multiplied by the number of parameters. Past generations are calculated by deviding by the number of chains per population}
\item{useMultiT}{whether to downscale Temperature of result components during burnin}
}

}}

\describe{ \item{Initial state: \code{X}}{
If initial state X is not specified, the last column (generation) of Z is utilized.
If in addition to X, logDenX is specified as a numeric vector, the fLogDen will not be avaluated for the initial state of the chains.
All the results of fLogDen for the initial state must be finite.
}}

\describe{ \item{Acceptance rate}{
The acceptance rate is tracked for each chain across ctrl$pAcceptWindowWidth generations. \cr
If acceptance rates drops to low values, this might be because of bad localization
,i.e the requirement to use states from more distant past.
In order to improve localization, less parameters or more chains per population are required.
}}}
\value{a list with entries of populations, each entry is a list
\item{thin}{thinning interval that has been used}
\item{dInfos}{list of information on densities (argument \code{dInfos})}
\item{blocks}{list of information on blocks (argument \code{blocks})}
\item{YPos}{list of column positions in Y, a list with entries \describe{
\item{accepted}{integer vector of positions of acceptance indication of block at given index}
\item{resLogDen0}{integer scalar: postion before first column of results of fLogDen}
}}
\item{pops}{info on each population. A list with entries: \describe{
\item{upperParBounds}{upper parameter bounds for sampling}
\item{lowerParBounds}{lower parameter bounds for sampling}
\item{splits}{named numeric vector: splitting history}
\item{spaceInd}{the space replicate that the population belongs to}
\item{parms}{numeric array (steps x parms x chains): collected states, including the initial states}
\item{temp}{numeric vector (nSample+1): global temperature, i.e. cost reduction factor}
\item{pAccept}{acceptance rate of chains (nStep x nChainPop)}
\item{resLogDen}{numeric array (steps x resComps x chains): results components of fLogDen of blocks}
\item{logDen}{numberic array (steps x iDen x chains): results summed over blocks}
}}}

\author{Thomas Wutzler}



\seealso{\code{\link{calcDEMCTemp}}
\code{\link{logDenGaussian}}}
\examples{
data(twLinreg1); attach( twLinreg1 ) 

.nPop=2
Zinit <- ZinitPops <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChains=4*.nPop, nPop=.nPop)
pops <- list(
	pop1 <- list(
		Zinit = ZinitPops[1:3,,1:4,drop=FALSE]	# the first population with less initial conditions
		,nGen=10
	),
	pop2 <- list(
		Zinit = ZinitPops[,,5:8,drop=FALSE]	# the first population with less initial conditions
		,nGen=15
		,T0=10
	)
)
#tmp <- .checkPop(pops[[1]])

# both blocks compare the same model against the same observations
blockDefault <- list(
	fLogDen=logDenGaussian
	,resCompNames=c("obs","parms")
	,TFix=c(parms=1)
	,argsFLogDen = list(
		fModel=dummyTwDEMCModel,		### the model function, which predicts the output based on theta 
		obs=obs,				### vector of data to compare with
		invCovar=invCovar,		### the inverse of the Covariance of obs (its uncertainty)
		xval=xval
	)
)
# first block updates only parameter a, second block updates only parameter b
# contrary to this example, correlated parameters should not be distributed across blocks
blocks <- list(
	blockA = within(blockDefault,{
			compPos="a"
			argsFLogDen <- c( argsFLogDen, list(
					thetaPrior= thetaTrue["a"]	### the prior estimate of the parameters
					,invCovarTheta = invCovarTheta[1,1,drop=FALSE]	### the inverse of the Covariance of the prior parameter estimates
					,blockIndices=1
				)) 
		})
	,blockB = within(blockDefault,{
			compPos="b"
			argsFLogDen <- c( argsFLogDen, list(
					thetaPrior= thetaTrue["b"]	### the prior estimate of the parameters
					,invCovarTheta = invCovarTheta[2,2,drop=FALSE]	### the inverse of the Covariance of the prior parameter estimates
					,blockIndices=2
				)) 
		})
)
#tmp <- .checkBlock(blocks[[2]],parNames=c("a","b"))
block <- blocks[[2]] 
block <- blocks[[1]] 
#mtrace(logDenGaussian)
(tmp <- do.call( logDenGaussian, c(list(pops[[1]]$Zinit[1,,1]), block$argsFLogDen) ))


#mtrace(twDEMCBlockInt)
#mtrace(updateBlockTwDEMC)
res <- twDEMCBlockInt( pops=pops, blocks=blocks, nGen=60)
str(res[[2]])
#plot(res[[2]]$temp)

# load the restart file and continue
load( file=restartFilename )	# variable resRestart.twDEMC
getNGen(resRestart.twDEMC)		# 48, last thinnging interval (each thin=4) before 50
res2 <- twDEMCBatch(resRestart.twDEMC)	# continue without needing to respecify parameters
getNGen(res2)					# 60 as nGen
res3 <- twDEMCBatch(res2, nGen=100)	    # continue even further without needing to respecify parameters
getNGen(res3)					# 100 as nGen

detach()

}
