\name{popMeansTwDEMC}
\alias{popMeansTwDEMC}
\title{popMeansTwDEMC}
\description{Calculating population means across chains within each population, and smooth time series.}
\usage{popMeansTwDEMC(x, nPop, kSmooth = NULL)}
\arguments{
  \item{x}{a matrix with columns chains}
  \item{nPop}{number of populations}
  \item{kSmooth}{weights to the filter function, or just number of points}
}



\author{Thomas Wutzler}



\seealso{\code{\link{popApplyTwDEMC}}
\code{\link{subChains.twDEMC}}}
\examples{
data(twdemcEx1)
ex1c <- concatPops(twdemcEx1)
# mean logDen for each case, i.e. step, by population
#mtrace(popMeansTwDEMC)
res1 <- popMeansTwDEMC( ex1c$logDen[,1,], nPop=getNPops(ex1c) )
matplot(res1)
# shifting mean across 4 cases
res2 <- popMeansTwDEMC( ex1c$logDen[,1,], nPop=getNPops(ex1c), kSmooth=5 )
matplot(res2, type="l", add=TRUE)
}
