\name{replaceZinitCases}
\alias{replaceZinitCases}
\title{replaceZinitCases}
\description{Replaces states of Zinit by sampling the other states that are marked good.}
\usage{
replaceZinitCases(Zinit, boMat)
}
\arguments{
  \item{Zinit}{initial states of the form required by \code{\link{twDEMCInt}}}
  \item{boMat}{boolean matrix specifying good cases with rows cases and columns chains. If it is a vector, then matrix is constructed number of chains last dimension}
}
\details{Samples for the first half of chains are sampled from good cases of the second half of chains.
Samples for the second half of chains are sampled from good cases of first half of chains.}
\value{Zinit, with several cols (parameter vectors) replaced by other cols}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}
\code{\link{replaceZinitNonFiniteLogDens}}}

