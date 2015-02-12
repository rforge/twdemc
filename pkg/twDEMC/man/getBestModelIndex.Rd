\name{getBestModelIndex}
\alias{getBestModelIndex}
\title{getBestModelIndex}
\description{select the best model based on (temperated) logDensity components}
\usage{getBestModelIndex(logDenT, dInfos)}\arguments{
  \item{logDenT}{numeric matrix (nStep x nResComp): logDensity (highest are best)}
  \item{dInfos}{list of lists with entry resCompPos (integer vector) specifying the position of result components for each density}
}



\author{Thomas Wutzler}



\seealso{\code{\link{orderLogDen}}}
\examples{
logDenT <- cbind( -sample(5)/2, -sample(5), -sample(5) )
#dInfos <- list( d1=list(resCompPos=1:2), d2=list(resCompPos=3) )
dInfos <- list( d1=list(resCompPos=2), d2=list(resCompPos=3) )
getBestModelIndex(logDenT, dInfos)
-logDenT
}
