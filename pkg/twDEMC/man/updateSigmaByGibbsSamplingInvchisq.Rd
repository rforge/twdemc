\name{updateSigmaByGibbsSamplingInvchisq}
\alias{updateSigmaByGibbsSamplingInvchisq}
\title{updateSigmaByGibbsSamplingInvchisq}
\description{update Variance by sampling from a scaled inverse Chi-square distribution with using prior information on sigma }
\usage{updateSigmaByGibbsSamplingInvchisq(theta, argsFUpdateBlock = NULL, 
    fResid, nu0 = 0, sigma20, ...)}\arguments{
  \item{theta}{numeric vector: current state that is used in density function of the block}
  \item{argsFUpdateBlock}{list with possible component \code{intermediate} that will be provided to \code{fResid}}
  \item{fResid}{function to calculate residuals between prediction and observations, must return a numeric vector}
  \item{nu0}{weight of prior: equivalent number of observations, defaults to no prior (zero obs)}
  \item{sigma20}{prior estimate of variance, needs to be specified only when giving  nu0 different from zero}
  \item{\dots}{further arguemnts to fResid}
}
\details{Gelman et al., 2003, p.51: Estimating variance with known mean.

See the vignette on unknown observation uncertainty.}
\value{list with components
\item{accepted}{boolean scalar: if step was accepted}
\item{xC}{numeric vector: components of position in parameter space that are being updated}
\item{intermediate}{new intermediate state, if provided by fResid
}}

\author{Thomas Wutzler}





