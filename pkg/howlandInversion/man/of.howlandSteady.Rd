\name{of.howlandSteady}
\alias{of.howlandSteady}
\title{of howlandSteady}
\description{Objective function for comparing against Howland data, assuming constant input and steady state C-stocks, which determines k-Values.}
\usage{
of.howlandSteady(normpopt, logLikAccept = numeric(0), metropolisStepTemp = structure(rep(1, length(logLikAccept)), names = names(logLikAccept)), ..., model, poptDistr, obs = Howland14C$obsNutrientSite, times = 1950:2007, input = Howland14C$litter, parms = HowlandParameterPriors$parms0, fCalcBiasedObs = NULL, argsFCalcBiasedObs = list(), fCalcBiasedInput = meanInput, argsFCalcBiasedInput = list(), popt.names = names(normpopt), fTransOrigPopt = transOrigPopt.default, doStopOnError = FALSE, includeStreams = c(names(obs), "parms"), delta14Catm = c14Constants$delta14Catm, fCalcIROLayer = calcIROLayer)
}
\arguments{
  \item{normpopt}{numeric vector: the point in normal parameter space}
  \item{logLikAccept}{numeric named vector the logLikelihood of previously accepted parameters, see details}
  \item{metropolisStepTemp}{the temperature for current metropolis decision}
  \item{\dots}{further parameters to model$fSolve}
  \item{model}{model to evaluate. A list with following components \describe{
\item{modMeta}{meta information about the model. see \code{\link{twCreateModMeta}} }
\item{fInitState}{function to initialize state based on parameters and amendment. see \code{\link{initStateSoilMod}} }
\item{fSolve}{function to ODE across time. see \code{\link{solveICBM1}} }
}}
  \item{poptDistr}{information on parameter distributions, list with entries \describe{
\item{trans}{character vector: type of distribtution (norm,lognorm,logitnorm)}
\item{mu}{numeric vector: distribution parameter mu, i.e. expected values at normal scale}
\item{invsigma}{numeric matrix: inverse of the distribution parameter sigma on multivariate (transformed) normal distribution}
}
See also \code{\link{twQuantiles2Coef}}.
Parameters order must correspond to normpopt. Use \code{\link{twConstrainPoptDistr}} to select a constrained subset from all parameters.}
  \item{obs}{see \code{\link{Howland14C}}}
  \item{times}{time points in yr to be modelled and compared.}
  \item{input}{carbon inputs to reservoir C see \code{\link{Howland14C}$litter}.}
  \item{parms}{default parameters (for the non-optimized ones)}
  \item{fCalcBiasedObs}{function(obs,padj,...){obs} possibility to account for bias and to optimize bias parameters}
  \item{argsFCalcBiasedObs}{further arguments to fCalcBiasedObs}
  \item{fCalcBiasedInput}{function(input,padj,...){obs} possibility to account for bias and to optimize bias parameters}
  \item{argsFCalcBiasedInput}{further arguments to fCalcBiasedInput}
  \item{popt.names}{names of the parameters. They are sometimes stripped by fitting algorithms.}
  \item{fTransOrigPopt}{function that translates parameters from normal to original scale}
  \item{doStopOnError}{by default -Inf is returned on error, set to TRUE for debugging}
  \item{includeStreams}{character vector of subset of data streams to include}
  \item{delta14Catm}{14C signature of the atmosphere: signature of litter inputs after a time lag.}
  \item{fCalcIROLayer}{function to calculate iRofO-Layer}
}
\details{objective function for fitting SoilMod_FS to respRate timeseries of Hamer incubation experiment of both control, amendm and c14obs
inclDataSeries can be set to a subset of 1:3 corresponding to control, amendment_large, and amendment_small respectively to save computing time

Supports a multi-step Metropolis descision. If \code{logLikAccept["parms"]} is provided, 
then a Metropolis descision is done based only on the parameters.
If it fails, then -Inf is returned. 
The possible costly evaluation of fModel is avoided.
\cr
Similarly, the small amendmend experiment and the control are not calculated, if 
observations are rejected be the large amendment experiment.
Provide a subset of \code{logLikAccept[c("amdend","c14obs","amdend_small","c14obs_small")]} to support early rejection.}


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





