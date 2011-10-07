\name{calcDEMCTempProp}
\alias{calcDEMCTempProp}
\title{calcDEMCTempProp}
\description{Calculate Temperature of components }
\usage{
calcDEMCTempProp(temp, diffLogDen, rFracMin = 1/4)
}
\arguments{
  \item{temp}{the maximum temperature}
  \item{diffLogDen}{expected difference in LogDensitys expected-accepted per datastream}
  \item{rFracMin}{fraction of max DiffDensity  below which temperatue is scaled down to yield  larger importance}
}

\value{vector of temperatures corresponding to diffLogDen with maximum corresponding to temp}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





