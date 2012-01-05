\name{ggplotDensity.twDEMC}
\alias{ggplotDensity.twDEMC}
\title{ggplotDensity twDEMC}
\description{Plotting the densities for each parameter.}
\usage{
ggplotDensity.twDEMC(res,
poptDistr = NULL,
pMin = 0.05,
doTransOrig = FALSE,
doDispLogDen = TRUE)
}
\arguments{
  \item{res}{the twDEMC whose densities to plot}
  \item{poptDistr}{parameter Distributions for the prior, usually \code{poptDistr <- \link{twConstrainPoptDistr}(poptNames,HamerParameterPriors$parDistr )}}
  \item{pMin}{if > 0, the results are constrained to quantiles of rLogDen>percMin. Can avoid extremes}
  \item{doTransOrig}{if TRUE, parameters are translated to original scale}
  \item{doDispLogDen}{include density of LogDensitys}
}



\author{Thomas Wutzler}



\seealso{\code{\link{plotMarginal2D}}
\code{\link{twDEMCInt}}}

