\name{stackChainsPop.twDEMC}
\alias{stackChainsPop.twDEMC}
\title{stackChainsPop twDEMC}
\description{Combine MarkovChains of each population of a twDEMC. }
\usage{\method{stackChainsPop}{twDEMC}(x, ...)}\arguments{
  \item{x}{
}
  \item{\dots}{
}
}



\author{Thomas Wutzler}



\seealso{\code{\link{stackChains.twDEMC}}
\code{\link{subChains.twDEMC}}}
\examples{
data(twdemcEx1)
ex1c <- concatPops(twdemcEx1)
res <- stackChainsPop(ex1c)
str(res)
(tmp1 <- head(res[,,1]))
}
