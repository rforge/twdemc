\name{twDEMCBlockInt}
\alias{twDEMCBlockInt}
\title{twDEMCBlockInt}
\description{Differential Evolution Markov Chain with blocked parameter update}
\usage{twDEMCBlockInt(pops = list(list(parms = list(), nGen = 12, TProp = NULL, 
    X = NULL, logDenCompX = NULL, spaceInd = 1, upperParBounds = numeric(0), 
    lowerParBounds = numeric(0), splits = numeric(0))), dInfos = list(list(fLogDen = NULL, 
    xC = c(), logDenCompC = c(), parUpdateDenC = c(), compPosDen = 1:nrow(pops[[1]]$parms), 
    argsFLogDen = list(), intResComp = vector("integer", 0), 
    intermediate = "modOutput", nObs = c())), blocks = list(list(dInfoPos = 1, 
    compPos = dInfos[[1]]$compPosDen, fUpdateBlock = updateBlockTwDEMC, 
    argsFUpdate = list(), requiresUpdatedDen = TRUE)), TSpec = matrix(1, 
    ncol = 2, nrow = 1, dimnames = list(NULL, c("T0", "TEnd"))), 
    m0 = numeric(0), nGen = integer(0), spacePop = 1:nPop, controlTwDEMC = list(), 
    debugSequential = FALSE, remoteDumpfileBasename = NULL, fDiscrProp = NULL, 
    argsFDiscrProp = list(), doRecordProposals = FALSE, progressOutput = ".")}
\arguments{
  \item{pops}{list of population infos for each population, each info is a list with components\describe{
\item{parms}{list of matrices (nState x nParm  x nChain) initial states for each population see details and \code{\link{initZtwDEMCNormal}}.}
\item{nGen}{scalar integer: number of generations, if given overwrites \code{pops[[i]]$nGen}}
\item{TProp}{numeric vector (nResComp) temperature proportions of result components.
It can also be given as character vector with names of result components, however, be aware that htis fails if several logDen may return components of the same name}
\item{X}{numeric matrix (nParm x nChainPop) initial state}
\item{logDenCompX}{numeric matrix (nComp x nChain): logDen components of initial state X, see details}
\item{spaceInd}{the space replicate that this population belongs to}
\item{upperParBounds}{named numeric vectors: giving upper parameter bounds: lowerBound < par <= upperBound.
For exploring subspaces of the limiting distribution, see details}
\item{lowerParBounds}{similar to upperParBounds: sample > bound}
\item{splits}{named numeric vector of splitting points, used to remerge divided subspaces}
}}
  \item{dInfos}{named list of used density functions. Each entry is a list with components\describe{
\item{fLogDen}{\code{function(theta, ...)} calculates a vector of logDensities
corresponding to different data streams of parameter vector theta 
\cr or \code{function(theta, logDenCompX, metropolisStepTemp, ...)}}
\item{xC}{numeric vector: current accepted state}
\item{logDenCompC}{numeric vector: result components of fLogDen for current position}
\item{parUpdateDenC}{integer vector: logDensity that recently updated parameter at given index
to handle first steps in Multi-step Metropolis decision internally. 
See details.}
\item{compPosDen}{index or names of the parameter components that are used in this density function}
\item{argsFLogDen}{further arguments passed to fLogDen}
\item{intResComp}{integer or character vector: indices or names of results components of fLogDen
that are used for internal Metropolis decisions}
\item{intermediate}{string: identifier in list of intermediate results that are shared between densities. See details.}
\item{nObs}{integer vector (nRespComp): number of observations for each result component, used in control of overfitting}
}

list describing the logDensities, see \code{\link{twDEMCBlockInt}}}
  \item{blocks}{list of parameter blocks, each block is a list with entries\describe{
\item{dInfoPos}{name or position to \code{fLogDenInfo}. Several blocks may share the same density but update different parameters}
\item{compPos}{names or index of the parameter components to be updated}
\item{fUpdateBlock}{function to update the parameters.
\cr It must return a list with first three components xC, logDenCompC, and parUpdateDenC 
as described in \code{\link{updateBlockTwDEMC}}}
\item{argsFUpdate}{further arguments passed to fUpdate}
\item{requiresUpdatedDen}{if fUpdateBlock does not depend on current density result, then set this to FALSE and save some calculation time}
}}
  \item{TSpec}{numeric matrix (nResComp x 2): specifying Initial and End Temperature of each fLogDen result component.
If only one row is specified, the Temperature is taken for all result components}
  \item{m0}{scalar integer: number of samples in initial population, if length==0 will be calculated from number of chains and number of parameters}
  \item{nGen}{scalar integer: number of generations, if given overwrites \code{pops[[i]]$nGen}}
  \item{spacePop}{the space replicate that each population belongs to. By default assume only one population per space, overridden by entry in pops}
  \item{controlTwDEMC}{control parameters influencing the update and the convergens of the chains (see details)}
  \item{debugSequential}{if TRUE apply is used instead of sfApply, for easier debugging}
  \item{remoteDumpfileBasename}{the basename of a dumpfile that is created on error on remote process (see example section)}
  \item{fDiscrProp}{function applied to proposal, e.g. to round proposals to to discrete possible values function(theta,...)}
  \item{argsFDiscrProp}{further arguments to fDiscrProp}
  \item{doRecordProposals}{if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y}
  \item{progressOutput}{output after each thinning interval}
}
\details{This is the central method for applying a Differential Evolution Markov Chain given a function of 
an unnormalized probability density.
It is invoked usually by (\code{\link{twDEMCBlock.array}} or \code{\link{twDEMCBlock.twDEMCPops}})

\describe{ \item{recognized control parameters in \code{controlTwDEMC}: }{\describe{
\item{thin}{thinning interval, Default: 4}
\item{F}{related to multiplicative error (F2=F/sqrt(2*Npar), see eps.mult, Default: 2.38}
\item{pSnooker}{probability of a snooker update (others parallel updates), Default 0.1}
\item{pGamma1}{probability of jumping to state of another chain (different modes), Default 0.1}
\item{epsMult}{>0 gives d-dimensional variation around gamma.
It adds scaled uncorrelated noise to the proposal. Default: 0.2 \cr 
Its advantage over eps.add is 
that its effect scales with the differences of vectors in the population whereas eps.add does not. 
if the variance of a dimensions is close to 0, eps.mult gives smaller changes. \cr 
A uniformly distributed error, i.e. F2*runif(1+-epsMult*prop) multiplied to difference vector 
from parallel update }
\item{epsAdd}{>0 is needed to ensure that all positions in the space can be reached.
For targets without gaps, it can set small or even to 0. Default 0 \cr 
sd of normally distributed error added to proposal by parallel or snooker update. }
\item{pAcceptWindowWidth}{number of thinning intervals back over which the acceptance rate is calculated Default 8}
\item{probUpDir}{probability of direction between two states of increasing Density, Default 0.5
\cr Increasing this during burin may accelerate convergence}
\item{initialAcceptanceRate}{numeric matrix (nBlock x nChains) initially assumed acceptance rate. Default 0.25
\cr Used to calculate the number of generations backwards to sample from}
\item{DRgamma}{factor for reducing step length [0..1) in delayed rejection step, 0 means no DR step, Default 0}
\item{minPCompAcceptTempDecr}{if acceptance rate drops below minPCompAcceptTempDecr+0.02 times this level, employ delayed rejection (DR), Default 0.15}
\item{pIndStep}{assume state to be independent, after on average about those number of accepted steps, Default 1.5}
\item{nPastGen}{factor for determining the number of recent past states to sample during burnin. Default 10
It is multiplied by the number of parameters. Past generations are calculated 
by deviding by the number of chains per population}
\item{returnIntermediates}{set to FALSE if intermediate result is large and transferring it between slaves slows down calculation}
\item{useConditionalProposal}{if set to TRUE, poposal steps are generated conditional on the state of the other blocks}
\item{controlOverfittingMinNObs}{scalar positive integer: minimum number of observations where bias is controlled. Need at least about 20 observations to work.}
}

}}

\describe{ \item{Initial state: \code{X}}{
If initial state X is not specified, the last column (generation) of Z is utilized.
If in addition to X, logDenX is specified as a numeric vector, the fLogDen will not be avaluated for the initial state of the chains.
All the results of fLogDen for the initial state must be finite.
}}

\describe{ \item{Acceptance rate}{
The acceptance rate is tracked for each chain across ctrl$pAcceptWindowWidth thinning intervals. \cr
If acceptance rates drops to low values, this might be because of bad localization
,i.e the requirement to use states from more distant past.
In order to improve localization, less parameters or more chains per population are required.
}}

\describe{ \item{Overfitting control}{
Cost, i.e. -2*logDensity, below the number of observations indicates an overfitting to a data stream.

By supplying the number of observations for one result component with \code{dInfo$nObs} one can 
activate the control for this overfitting. This works only, however, if there are enough observations. 
The minimum numbe of observations where overfitting control is activated is given with 
\code{ctrl$controlOverfittingMinNObs}. By setting this to Infinity, no overfitting control is applied.

With overfitting control, there is a lower bound for the cost of one result compoent 
, i.e. an upper bound of the corresponding \code{dInfo$fLogDen} result component. Values above this   
bound are decreased to this value. 
}}

\describe{ \item{intermediate results}{
Different density functions can share intermediate results, so that these do not need to be recomputed.
For example two densities can be based on the same model predictions for given parameters, but compare different
subsets of the predictions against different observations. If for a given parameter vector, the model output has been
computed in one density, it does not need to be recomputed in the second density function.

This is an advanced option. Care must be taken that the intermediate is really the same between
densities.

Description of both densities should specify the same identifier string in argument \code{dInfo$intermediate}. 
The densities should return the model output in attribute "intermediate" with the result vector of logDensity
components. \code{twDEMC} ensures that this result is provided with argument \code{intermediate}
in another call to the density function. 
If the parameters have changed, an empty list will be given with this argument.

For an example see test case \code{ofMultiIntermediate} in file unitTests/runittwDEMC.R.
Using densities \code{\link{denSparsePrior}} and \code{\link{denRichPrior}}. 

When using parallel execution on multiple cores, the transfer of very big intermediates can cause
performance issues. Use argument \code{controlTwDEMC$returnIntermediate=FALSE} to ensure that 
intemediates are stored only within process boundaries. When executing the next interval or the next round of
updates of all blocks, the intermediate will then be the empty list initially.
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
\item{temp}{numeric array (nSample+1, nResComp): temperature, i.e. cost reduction factor in each step for each datastream}
\item{pAccept}{acceptance rate of chains (nStep x nChainPop)}
\item{resLogDen}{numeric array (steps x resComps x chains): results components of fLogDen of blocks}
\item{logDen}{numberic array (steps x iDen x chains): results summed over blocks}
}}}

\author{Thomas Wutzler}



\seealso{\code{\link{calcDEMCTemp}}
\code{\link{logDenGaussian}}}
\examples{
    # see unit test 
    # twUtestF("twDEMCPops","distinctLogDen")
data(twLinreg1); attach( twLinreg1 ) 

# tracing error in remote session:
# Provide argument \code{remoteDumpfileBasename="dumpFile"}
# Then a dumpfile is created on remote error by \code{\link{sfRemoteWrapper}}.
# In order to trace the density function, the following can be done
if( FALSE ){
	.remoteDumpfileBasename="dumpfile"
	.remoteDumpfile <- paste(.remoteDumpfileBasename,".rda",sep="")
	load(.remoteDumpfile)
	debugger(get(.remoteDumpfileBasename))
	# choose last step (18)
	#require(debug)
	fDen <- remoteFunArgs$argsUpdateBlocksTwDEMC$argsFUpdateBlocks[[1]]$fLogDen
	mtrace(fDen); remoteFunArgs$argsUpdateBlocksTwDEMC$argsFUpdateBlocks[[1]]$fLogDen <- fDen
	do.call(remoteFun, c(remoteFunArgs, list(...)))	
	# go()
	# qqq()
}

}
