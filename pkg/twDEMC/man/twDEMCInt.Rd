\name{twDEMCInt}
\alias{twDEMCInt}
\alias{twDEMC}
\title{twDEMCInt}
\description{Differential Evolution Markov Chain}
\usage{
twDEMCInt(Zinit, nGen = 10, ..., controlTwDEMC = list(), X = NULL, logLikX = NULL, fLogLik, argsFLogLik = list(), resFLogLikX = NULL, intResCompNames = character(0), fLogLikScale = 1, debugSequential = FALSE, remoteDumpfileBasename = NULL, nPops = 1, fDiscrProp = NULL, argsFDiscrProp = list(), doRecordProposals = FALSE)
}
\arguments{
  \item{Zinit}{initial population: a matrix of number of parameters by number of individuals (d x m0 x nChains) see details and \code{\link{initZtwDEMCNormal}}.}
  \item{nGen}{number of generations, i.e. steps of the chain}
  \item{\dots}{further arguments passed to fLogLik}
  \item{controlTwDEMC}{DEMCzsControl parameters influencing the update and the convergens of the chains (see details)	}
  \item{X}{initial active population that will be evolved: a matrix of number of parameters by number of chains (d x N)
  if null, it will be initialized to the last row of Zinit}
  \item{logLikX}{log-Likelihood of initial values, may save time if already calculated}
  \item{fLogLik}{\code{function(theta, ...)} calculates a vector of log-Likelihoods corresponding to different data streams of parameter vector theta 
or \code{function(theta, resFLogLikX, metropolisStepTemp, ...)} to handle first steps in Multi-step Metropolis decision internally. See details.  }
  \item{argsFLogLik}{further arguments passed to fLogLik}
  \item{resFLogLikX}{numeric matrix: logLik components x chains, see details}
  \item{intResCompNames}{character vector: names of results components of fLogLik that are used for internal Metropolis decisions }
  \item{fLogLikScale}{scalar multiplied to the result of fLogLik 
  allows using functions of negative LogLikelihood (-1) or Gaussian misfit function (-1/2) instead of logLikelihood}
  \item{debugSequential}{if TRUE apply is used instead of sfApply, for easier debugging}
  \item{remoteDumpfileBasename}{the basename of a dumpfile that is created on error on remote process }
  \item{nPops}{the number of populations that do not mix: must be a factor of dimension nChains of Zinit}
  \item{fDiscrProp}{function applied to proposal, e.g. to round proposals to to discrete possible values function(theta,...)}
  \item{argsFDiscrProp}{further arguments to fDiscrProp}
  \item{doRecordProposals}{if TRUE then an array of each proposal together with the results of fLogLik are recorded and returned in component Y}
}
\details{This is the central method for applying a Differential Evolution Markov Chain.
It is invoked usually by \code{twDEMC} (\code{\link{twDEMC.array}} or \code{\link{twDEMCBatchInt}})

This method is based on Code of ter Braak  ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain 
with snooker updater and fewer chains. Statistics and Computing
http://dx.doi.org/10.1007/s11222-008-9104-9 .

\describe{ \item{Initial population: \code{Z}}{  
a matrix of number of parameters by number of individuals (d x m0 x Npop) \describe{
\item{d}{number of dimensions of parameter vector theta}
\item{m0}{initials : default 10d/Npop (smaller for fewer parameters) all drawn from prior}
\item{nChains}{number of chains}}
alternatively Zinit may be an twDEMC object, then it will be extended by nGen with parameters of former call
}}

\describe{ \item{Several populations: \code{nPops}}{
Chains within population are not independent.
In order to assess convergence, one must run several independent populations.
This is supported by specifying argument nPops. Then the n chains, i.e. dim(Zinit)[3], are grouped into several populations.
}}

\describe{ \item{Detailed control parameters: \code{controlTwDEMC}}{\describe{
\item{F}{related to multiplicative error (F2=F/sqrt(2*Npar), see eps.mult}
\item{pSnooker}{probability of a snooker update (others parallel updates)}
\item{pGamma1}{probability of jumping to state of another chain (different modes)}
\item{epsMult}{>0 gives d-dimensional variation around gamma. It adds scaled uncorrelated noise to the proposal. Its advantage over eps.add is that its effect scales with the differences of vectors in the population whereas eps.add does not. if the variance of a dimensions is close to 0, eps.mult gives smaller changes. \cr A uniformly distributed error, i.e. F2*runif(1+-epsMult*prop) multiplied to difference vector from parallel update}
\item{epsAdd}{>0 is needed to ensure that all positions in the space can be reached. For targets without gaps, it can set small or even to 0. \cr sd of normally distributed error added to proposal by parallel or snooker update.}
\item{thin}{thinning interval}
\item{T0}{initial and end Temperature to flatten the likelihood surface, for each population}
\item{useMultiT}{if TRUE Temperature for different data streams are scaled by logLik of accepted state}
\item{TFix}{named numeric vector of Temperatures of fLogLik components, whose Temperature is fixed}
\item{pAcceptWindowWidth}{number of generations back over which the acceptance rate is calculated}
\item{probUpDir}{probability of direction between two states of increasing Likelihood (increase during burin may accelerate convergence)}
\item{initialAcceptanceRate}{numeric vector (nPops) initially assumed acceptance rate. Used to calculate the number of generations backwards to sample from
}}}
}

\describe{ \item{Multistep Metropolis decisions: \code{intResCompNames} }{
For performance reasons one may decide to do a Metropolis decision on parameters
before evaluating the model in order to save the possibly costly model evaluation.
\cr
Argument \code{intResCompNames} specifies a character vector of return components that are handled internally in the objective function.
The first currently accepted value provided to next 
evaluation of the objective function is taken from \code{resFLogLikX}. If this 
argument is NULL, then a proper matrix of \code{resFLogLikX} will be initialized to -Inf, i.e. accepting the next step
\cr 
See the code of \code{\link{logLikGaussian}} for an example multistep version of fLogLik
}}

\describe{ \item{Initial state: \code{X}}{
If initial state X is not specified, the last column (generation) of Z is utilized.
If in addition to X, logLikX is specified, the fLogLik will not be avaluated for the initial state of the chains.
All the results of fLogLik for the initial state must be finite.
}}

\describe{ \item{Acceptance rate}{
The acceptance rate is tracked for each chain across ctrl$pAcceptWindowWidth generations.
\cr If acceptance rates drops to low values, this might be because of bad localization
,i.e the requirement to use states from more distant past.
In order to improve localization, less parameters or more chains per population are required.
}}}
\value{list of class \code{twDEMC} (with \code{nStep = M0+nGen\%/\%thin}) \describe{
\item{parms}{ array (d x nStep x nChain) the initial parameters (last row of M0) and the accepted parameter combinations}
\item{rLogLik}{ array (nStep x nChain) the logLik of the accepted parameter combinations}
\item{pAccept}{ array (nStep x nChain) the acceptance probability over previous ctrl$pAcceptWindowWidth steps }
\item{temp}{ vector (nStep) the Temperature used at the step }
\item{resFLogLikX}{ array (nInternalResult x nChain) of fLogLik results of currently accepted step}
\item{thin}{ numeric: thinning interval }
\item{Y}{ numeric matrix: record of all proposals together with results components and columns "rLogLik" and boolean "accepted".
	The length is the thinning intervals of the last 128steps or (if doRecordProposals) all the steps +1 for the first line of the initial state.}
}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{Further functionality of the twDEMC package deals with \itemize{
\item{ Batch invokation of twDEMC and storing intermediate results: \code{\link{twDEMCBatchInt}}  } 
\item{ Generating an initial population for twDEMC: \code{\link{initZtwDEMCNormal}}  }
\item{ Transforming the results of twDEMC: \code{\link{subChains.twDEMC}}  }
\item{ Transforming the parameter space: \code{\link{transOrigPopt.default}}  }
\item{ Invoking fLoglik with proposal in a parallel load balanced way: \code{\link{twCalcLogLikPar}}  }
\item{ Plotting routines: \code{\link{plotMarginal2D}}  }
\item{ Calculating marginal aggregates: \code{\link{marginals1d}}  }
\item{ Diagnostics for MCMC samples: \code{\link{checkConvergenceGelman}}  }
}
\code{\link{calcDEMCTemp}}
\code{\link{logLikGaussian}}
\code{\link{.generateXPropThin}}
\code{\link{.doDEMCSteps}}
\code{\link{.doDEMCStep}}}
\examples{
data(twLinreg1); attach(twLinreg1); plot(obs~xval); abline(theta0)

# setup the model and the function of Log-Likelihood of parameters
# ?dummyTwDEMCModel
# ?logLikGaussian
argsFLogLik <- list( fModel=dummyTwDEMCModel,obs=obs,invCovar=invCovar,	xval=xval )
do.call( logLikGaussian, c(list(theta=theta0),argsFLogLik)) #test calling the logLik function

# create an initial distribution of states around the estimate of a usual regression 
lmDummy <- lm( obs ~ xval, weights=1/sdObs^2)		# results without priors
(.expTheta <- structure(coef(lmDummy),names=c("a","b")) )
(.expCovTheta <- {tmp<-vcov(lmDummy); dimnames(tmp)<-list(names(.expTheta),names(.expTheta));tmp} )		# a is very weak constrained, negative covariance between a an b
confint(lmDummy)
.nPops=2
Zinit <- initZtwDEMCNormal( .expTheta, .expCovTheta, nChains=4*.nPops, nPops=.nPops)

# run the chains
res <-  twDEMC( Zinit, nGen=500, fLogLik=logLikGaussian, argsFLogLik=argsFLogLik, nPops=.nPops )

# plot the results
plotThinned(as.mcmc.list(res))	# no apparent burnin

# get the empirical standard deviation and 95\% confidence intervals from the DEMC-sample
sample <- stackChains(res)
apply(sample,2, sd)	
apply(sample,2, quantile, probs=c(0.025,0.975) )

# plot likelihood surface 
pairs(sample)
ds <- as.data.frame(apply(sample[ sample[,"rLogLik"] >= max(sample[,"rLogLik"]-1.9), ],2,function(var){
		grain <- diff(range(var))/60
		round(var/grain)*grain
	}))
tmpx <- sort(unique(ds$a)); tmpy <- sort(unique(ds$b))
dsog <- expand.grid(x=tmpx[-1], y=tmpy[-1])
dsog$z <- apply(as.matrix(dsog[,1:2]),1,function(xy){ 
		dss <- subset(ds,a==xy[1] & b==xy[2])
		if( 0<nrow(dss)) max(dss$rLogLik) else NA
	})
image( tmpx, tmpy,  matrix(dsog$z,nrow=length(tmpx)-1), col = rev(heat.colors(100)), xlab="a", ylab="b" )

#nicer and faster with packages lattice and Rcmdr
.tmp.f <- function(){
	library(lattice)
	# round numbers to see something in levelplot else points get too small
	sampleSig <- apply(sample[ sample[,"rLogLik"] >= max(sample[,"rLogLik"]-1.9), ],2,function(var){
			grain <- diff(range(var))/150
			round(var/grain)*grain
		})
	levelplot(rLogLik~a*b, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	
	library(Rcmdr)
	sampleSig <- sample[ sample[,"rLogLik"] >= max(sample[,"rLogLik"]-1.9), ]
	ds <- as.data.frame(sampleSig)
	scatter3d(ds$a, ds$rLogLik, ds$b
			, surface=FALSE
		   ,bg="white", axis.scales=TRUE, grid=TRUE, ellipsoid=FALSE, xlab="a" 
		   ,ylab="rLogLik", zlab="b"
	   	   , point.col=rev(heat.colors(100))[round(rescale(ds$rLogLik,to=c(1,100)))]
	)
}

# for mor examples with prior and simulated annealing see the vignettes
}
