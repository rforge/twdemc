\name{calcDEMCTempProp}
\alias{calcDEMCTempProp}
\title{calcDEMCTempProp}
\description{Calculate Temperature of components }
\usage{
calcDEMCTempProp(temp, diffLogLik, rFracMin = 1/4)
}
\arguments{
  \item{temp}{the maximum temperature}
  \item{diffLogLik}{expected difference in Log-Likelihoods expected-accepted per datastream}
  \item{rFracMin}{fraction of max DiffLikelihood  below which temperatue is scaled down to yield  larger importance}
}

\value{vector of temperatures corresponding to diffLogLik with maximum corresponding to temp}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





