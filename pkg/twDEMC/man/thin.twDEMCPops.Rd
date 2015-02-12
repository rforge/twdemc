\name{thin.twDEMCPops}
\alias{thin.twDEMCPops}
\title{thin twDEMCPops}
\description{Reduces the rows of an twDEMCPops object (list returned by \code{\link{twDEMCBlockInt}}) to correspond to a thinning of \code{newThin}.}
\usage{\method{thin}{twDEMCPops}(x, newThin = x$thin, start = 0, end = NA, ..., 
    doKeepBatchCall = FALSE)}\arguments{
  \item{x}{the twDEMCPops list to thin}
  \item{newThin}{finite numeric scalar: the target thinning factor, must be positive multiple of x$thin}
  \item{start}{numeric scalar: the start time of the chain. 
Note that time starts from zero.
If a vector or matrix is supplied (e.g. nGenBurnin) then the maximum is used}
  \item{end}{numeric vector the maximum end time of the populations. 
Note that time starts from zero.
If a scalar is supplied, then it is repeateed NPop times }
  \item{\dots}{
}
  \item{doKeepBatchCall}{wheter to retain the batch call attribute of x}
}



\author{Thomas Wutzler}



\seealso{\code{\link{subset.twDEMCPops}}}
\examples{
if( FALSE ){
	#mtrace(thin.twDEMCPops)
	tmp <- thin(res, start=4, newThin=8 )
	all.equals( c(1, 11), getNSamples(tmp))
}
data(twdemcEx1)
x <- twdemcEx1
c( nGen=getNGen(twdemcEx1), thin=twdemcEx1$thin, nSample=getNSamples(twdemcEx1) )

    .nGenBurnin <- max(getNGen(twdemcEx1))-50
thinned <- thin(twdemcEx1, start=.nGenBurnin)	# removing burnin period
c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned) )	#13 sample describing 48 generations

thinned <- thin(twdemcEx1, start=.nGenBurnin, newThin=8)	
c( nGen=getNGen(thinned), thin=thinned$thin, nSample=getNSamples(thinned), nGenBurnin=thinned$nGenBurnin )	#7 samples describing 48 generations
}
