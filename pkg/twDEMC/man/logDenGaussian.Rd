\name{logDenGaussian}
\alias{logDenGaussian}
\title{logDenGaussian}
\description{Invokes the model and calculate an unnormalized logDensity (-1/2*misfit) assuming (multivariate) Gaussian errors in both data in priors. }
\usage{logDenGaussian(theta, logDenAccept = numeric(0), metropolisStepTemp = c(parms = 1), 
    ..., fModel, theta0 = theta, obs, invCovar, thetaPrior = NULL, 
    invCovarTheta = NULL, namesTheta = NULL, blockIndices = NULL, 
    scale = -1/2)}
\arguments{
  \item{theta}{the parameter vector to be optimized (may be more than updated, if used with blocks).}
  \item{logDenAccept}{scalar: logDen for parms from revious run for two step Metropolis decision}
  \item{metropolisStepTemp}{numeric named vector: the temperature for internal metropolis step}
  \item{\dots}{any other arguments passed to fModel}
  \item{fModel}{the model function, which predicts the output based on theta}
  \item{theta0}{parameter vector, first argument to fModel. Before invocation components theta overwrite theta0}
  \item{obs}{vector of data to compare with}
  \item{invCovar}{the inverse of the Covariance of obs (its uncertainty)
<< alternatively a vector of variances (diagonal covariance matrix) can be supplied and calculation is much more efficient}
  \item{thetaPrior}{the prior estimate of the parameters}
  \item{invCovarTheta}{the inverse of the Covariance of the prior parameter estimates
<< alternatively a vector of variances (diagonal covariance matrix) can be supplied and calculation is much more efficient}
  \item{namesTheta}{names assigned to theta (if not NULL), before invoking mofModel}
  \item{blockIndices}{integer vector: index of the components in theta and theta0 that should be regarded in this block}
  \item{scale}{factor to mulitply the misfit (e.g. -1/2 to obtain the unnormalized logDensity)}
}
\details{If thetaPrior is not specified (NULL) then no penalty is assigned to parameters.

Supports a two-step Metropolis descision. If \code{logDenAccept["parms"]} is provided, 
then a Metropolis descision is done based only on the parameters.
If it fails, then \code{c(obs=NA, parms=-Inf)} is returned. 
The possible costly evaluation of fModel is avoided.}
\value{the misfit: scale *( t(tmp.diffObs) \%*\% invCovar \%*\% tmp.diffObs + t(tmp.diffParms) \%*\% invCovarTheta \%*\% tmp.diffParms )}

\author{Thomas Wutzler}



\seealso{\code{\link{twDEMCBlockInt}}
\code{\link{dummyTwDEMCModel}}}

