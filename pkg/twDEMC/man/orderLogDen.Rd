\name{orderLogDen}
\alias{orderLogDen}
\title{orderLogDen}
\description{rank the log-Densities, starting with highest densities}
\usage{orderLogDen(ss, nDen = attr(ss, "nBlock"), decreasing = TRUE)}\arguments{
  \item{ss}{numeric matrix ( nRec x nDen+nPar )}
  \item{nDen}{number of densities}
  \item{decreasing}{argument to order, to start with the highest ranks, i.e best models}
}



\author{Thomas Wutzler}



\seealso{\code{\link{getBestModelIndices}}
\code{\link{sumLogDenCompBlocks}}}
\examples{
if( FALSE ){    # no twDEMCSA result available for now
    #assume resBlock is a result of twDEMCSA with multiple densities
    ss <- stackChains(resBlock, useTemperatedLogDen=TRUE)         # get the temperated summed densities for each block
    ss[ head(orderLogDen(ss)), ]     # display the best results
}
}
