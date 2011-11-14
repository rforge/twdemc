\name{ggplotDensity.poptDistr}
\alias{ggplotDensity.poptDistr}
\title{ggplotDensity poptDistr}
\description{Plotting the densities for each parameter.}
\usage{
ggplotDensity.poptDistr(poptDistr, pMin = 0.005, parmsBounds = NULL, plotUpperQuantile = TRUE, doTransOrig = TRUE)
}
\arguments{
  \item{poptDistr}{parameter Distributions for the prior, usually \code{poptDistr <- twConstrainPoptDistr(poptNames,HamerParameterPriors$parDistr )}}
  \item{pMin}{range of the distribution from pMin to 1-pMin}
  \item{parmsBounds}{list parName <- c(mode, upperQuantile)}
  \item{plotUpperQuantile}{wheter to include upper quantile (set to FALSE if this inflates the displayed domain)}
  \item{doTransOrig}{set to FALSE to display transform to normal scale}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{plotMarginal2D}}
\code{\link{twDEMCInt}}}

