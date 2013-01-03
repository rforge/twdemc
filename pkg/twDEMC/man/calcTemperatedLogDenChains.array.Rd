\name{calcTemperatedLogDenChains.array}
\alias{calcTemperatedLogDenChains.array}
\title{calcTemperatedLogDenChains array}
\description{Rescale Log-Densities by given Temperatures within chains}
\usage{\method{calcTemperatedLogDenChains}{array}(x, temp, ...)}
\arguments{
  \item{x}{numeric matrix (nStep x nResComp x nChain), e.g. pop$resLogDen}
  \item{temp}{numeric vector (nResComp): Temperature, i.e. cost reduction factor}
  \item{\dots}{
}
}

\value{numeric array (nStep x nComp x nChain): rescaled Log-Density for each chain.}

\author{Thomas Wutzler}



\seealso{\code{\link{sumLogDenComp}}}
\examples{FALSE}
