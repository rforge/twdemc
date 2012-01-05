\name{concatPops}
\alias{concatPops}
\title{concatPops}
\description{Generic method behind \code{\link{concatPops.twDEMCPops}}}
\usage{
concatPops(...)
}
\arguments{
  \item{\dots}{
}
}



\author{Thomas Wutzler}




\examples{
if( FALSE ){
	getNSamples(tmp <- concatPops(res))
	getNChains(tmp)
	getNPops(tmp)
	#mtrace(concatPops.twDEMCPops)
	getNSamples(tmp <- concatPops(res,minPopLength=10))
	getNChains(tmp)
	getNPops(tmp)
}
}
