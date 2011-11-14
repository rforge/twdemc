\name{popApplyTwDEMC}
\alias{popApplyTwDEMC}
\title{popApplyTwDEMC}
\description{Applying a function across all chains of one population for each case.}
\usage{
popApplyTwDEMC(x, nPops, FUN, ...)
}
\arguments{
  \item{x}{a matrix with columns chains or array with last dimension chain}
  \item{nPops}{number of populations}
  \item{FUN}{function to apply to population submatrix}
  \item{\dots}{further arguemtns to FUN}
}

\value{array with last dimenstion correponding to population}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{popMeansTwDEMC}}
\code{\link{subChains.twDEMC}}}
\examples{
data(twdemcEx1)
# mean rLogDen for each case, i.e. step, by population
nPops=getNPops(twdemcEx1)
popApplyTwDEMC( twdemcEx1$rLogDen, nPops=nPops, apply, 1, mean )	
# stack rLogDen for each population
popApplyTwDEMC( twdemcEx1$rLogDen, nPops=nPops, as.vector )
#stack param columns by population
(tmp <- popApplyTwDEMC( twdemcEx1$parms, nPops=nPops, function(x){ abind(twListArrDim(x),along=2) }))	
}
