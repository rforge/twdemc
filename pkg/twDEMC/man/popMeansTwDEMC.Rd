\name{popMeansTwDEMC}
\alias{popMeansTwDEMC}
\title{popMeansTwDEMC}
\description{Calculating population means across chains, and smooth time series.}
\usage{
popMeansTwDEMC(x, nPops, kSmooth = NULL)
}
\arguments{
  \item{x}{a matrix with columns chains}
  \item{nPops}{number of populations}
  \item{kSmooth}{weights to the filter function, or just number of points}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{popApplyTwDEMC}}
\code{\link{subChains.twDEMC}}}
\examples{
data(twdemcEx1)
# mean rLogDen for each case, i.e. step, by population
res1 <- popMeansTwDEMC( twdemcEx1$rLogDen, nPops=ncol(twdemcEx1$temp) )
matplot(res1)
# shifting mean across 4 cases
res2 <- popMeansTwDEMC( twdemcEx1$rLogDen, nPops=ncol(twdemcEx1$temp), kSmooth=5 )
matplot(res2, type="l", add=TRUE)
}
