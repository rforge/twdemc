\name{thin.twDEMC}
\alias{thin.twDEMC}
\title{thin twDEMC}
\description{Reduces the rows of an twDEMC object (list returned by \code{\link{twDEMCBlockInt}}) to correspond to a thinning of \code{newThin}.}
\usage{\method{thin}{twDEMC}(x, newThin = x$thin, start = 0, end = NA, ..., doKeepBatchCall = FALSE)}
\arguments{
  \item{x}{the twDEMC list to thin}
  \item{newThin}{finite numeric scalar: the target thinning factor, must be positive multiple of x$thin}
  \item{start}{numeric scalar: the start time of the chain. 
Note that time starts from zero.
If a vector or matrix is supplied (e.g. nGenBurnin) then the maximum is used}
  \item{end}{numeric scalar: the maximum end time of the chains. 
Note that time starts from zero.
If a vector or matrix is supplied (e.g. nGenBurnin) then the maximum is used}
  \item{\dots}{
}
  \item{doKeepBatchCall}{wheter to retain the batch call attribute of x}
}



\author{Thomas Wutzler}



\seealso{\code{\link{subChains.twDEMC}}}
\examples{
data(twdemcEx1)
x <- twdemcEx1
c( nGen=getNGen(twdemcEx1), thin=twdemcEx1$thin, nSample=getNSamples(twdemcEx1) )

    .nGenBurnin <- max(getNGen(twdemcEx1))-70
thinned <- thin(twdemcEx1, start=.nGenBurnin)	# removing burnin period
c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin )	#18 sample describing 68 generations

thinned <- thin(twdemcEx1, start=.nGenBurnin, newThin=8)	
c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin )	#9 samples describing 64 generations
}
