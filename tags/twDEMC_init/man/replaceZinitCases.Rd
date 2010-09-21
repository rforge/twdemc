\name{replaceZinitCases}
\alias{replaceZinitCases}
\title{replaceZinitCases}
\description{Replaces states of Zinit that yield non-finite rLogLik by sampling other states.}
\usage{
replaceZinitCases(Zinit, boMat)
}
\arguments{
  \item{Zinit}{initial states of the form required by \code{\link{twDEMCInt}}}
  \item{boMat}{boolean matrix with rows cases and columns parameters}
}
\details{Samples for the first half of chains are sampled from good cases of the second half of chains.
Samples for the second half of chains are sampled from good cases of first half of chains.}
\value{Zinit, with several rows replaced by other rows}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}
\code{\link{replaceZinitNonFiniteLogLiks}}}

