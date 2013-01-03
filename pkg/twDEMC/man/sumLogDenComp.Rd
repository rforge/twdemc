\name{sumLogDenComp}
\alias{sumLogDenComp}
\title{sumLogDenComp}
\description{sum LogDen components within density}
\usage{sumLogDenComp(resLogDen, dInfos = NULL)}
\arguments{
  \item{resLogDen}{numeric matrix (nStep x nResComp x nChain): of log-Densities}
  \item{dInfos}{list of densities each a list with entry resCompPos: integer vector, specifying the result components within density
<< If dInfos is of length 0, all components are summed and a matrix instead of a array is returned}
}

\value{numeric array (nStep x nDensities x nChain )}

\author{Thomas Wutzler}



\seealso{\code{\link{calcTemperatedLogDenChains.array}}}

