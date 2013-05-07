\name{subsetF.twDEMCPops}
\alias{subsetF.twDEMCPops}
\title{subsetF twDEMCPops}
\description{keeps only cases within population which evaluate to TRUE for a given function.}
\usage{\method{subsetF}{twDEMCPops}(x, fKeep, ...)}
\arguments{
  \item{x}{object of class twDEMCPops}
  \item{fKeep}{function(pop) returning an boolean matrix (nStep x nChain) of cases to keep
<< ,alternatively returning an integer matrix (niStep x nChain) with the indices to keep
<< ,alternatively returning an integer or boolean vector, that is applied to each chain}
  \item{\dots}{further arguments to fKeep}
}
\details{The samples are redistributed across chains}


\author{Thomas Wutzler}




\examples{
data(twdemcEx1)
range(concatPops(twdemcEx1)$parms[,"a",]) # spanning 9 to 11
#pop <- twdemcEx1$pops[[1]]
# note that the numer of samples across chains within one population is allowed to differ
fKeep <- function(pop){ tmp <- (pop$parms[,"a",] < 10) }
res <- subsetF(twdemcEx1, fKeep )
plot( as.mcmc.list(res), smooth=FALSE )
getNSamples(res)
}
