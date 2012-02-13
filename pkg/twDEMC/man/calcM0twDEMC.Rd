\name{calcM0twDEMC}
\alias{calcM0twDEMC}
\title{calcM0twDEMC}
\description{Calculate appropriate number of cases for initializing twDEMC.}
\usage{calcM0twDEMC(nPar, nChainPop)}
\arguments{
  \item{nPar}{the number of parameters to estimate}
  \item{nChainPop}{the number of chains per population}
}
\details{see terBraak 2006 and 2008}
\value{length of each chain so that each population is initialized with 8*nPar cases }

\author{Thomas Wutzler}



\seealso{\code{\link{initZtwDEMCNormal}}}

