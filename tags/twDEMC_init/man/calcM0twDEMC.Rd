\name{calcM0twDEMC}
\alias{calcM0twDEMC}
\title{calcM0twDEMC}
\description{Calculate appropriate number of cases for initializing twDEMC.}
\usage{
calcM0twDEMC(nPar, nPops, nChains)
}
\arguments{
  \item{nPar}{the number of parameters to estimate}
  \item{nPops}{the number of independent populations}
  \item{nChains}{the number of chains}
}
\details{see terBraak 2006 and 2008}
\value{length of each chain so that each population is initialized with 8*nPar cases }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}}

