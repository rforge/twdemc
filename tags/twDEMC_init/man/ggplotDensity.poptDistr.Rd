\name{ggplotDensity.poptDistr}
\alias{ggplotDensity.poptDistr}
\title{ggplotDensity poptDistr}
\description{Plotting the densities for each parameter.}
\usage{
ggplotDensity.poptDistr(poptDistr, pMin = 0.005, parmsBounds = NULL)
}
\arguments{
  \item{poptDistr}{parameter Distributions for the prior, usually \code{poptDistr <- twConstrainPoptDistr(poptNames,HamerParameterPriors$parDistr )}}
  \item{pMin}{range of the distribution from pMin to 1-pMin}
  \item{parmsBounds}{list parName <- c(mode, upperQuantile)}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





