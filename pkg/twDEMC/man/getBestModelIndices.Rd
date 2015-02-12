\name{getBestModelIndices}
\alias{getBestModelIndices}
\title{getBestModelIndices}
\description{select the best models based on (temperated) logDensity components}
\usage{getBestModelIndices(resLogDen, dInfos, prob = 0.1, doOrder = FALSE)}\arguments{
  \item{resLogDen}{numeric matrix (nStep x nResComp): logDensity (highest are best)}
  \item{dInfos}{list of lists with entry resCompPos (integer vector) specifying the position of result components for each density}
  \item{prob}{proportion of good indices to return}
  \item{doOrder}{if set to TRUE results will be ordered, yield a slight performance penalty}
}



\author{Thomas Wutzler}



\seealso{\code{\link{orderLogDen}}}
\examples{
logDenT <- cbind( -sample(5)/2, -sample(5), -sample(5) )
#dInfos <- list( d1=list(resCompPos=1:2), d2=list(resCompPos=3) )
dInfos <- list( d1=list(resCompPos=2), d2=list(resCompPos=3) )
getBestModelIndices(logDenT, dInfos, 0.2)
}
