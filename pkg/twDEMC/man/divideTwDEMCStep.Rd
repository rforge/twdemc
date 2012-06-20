\name{divideTwDEMCStep}
\alias{divideTwDEMCStep}
\title{divideTwDEMCStep}
\description{run twDEMCBlock on subspaces}
\usage{
divideTwDEMCStep(aTwDEMC,
qPop = numeric(0),
nGen = 512,
controlTwDEMC = list(),
m0 = calcM0twDEMC(getNParms(aTwDEMC),
getNChainsPop(aTwDEMC)),
nRowsMin,
nRowsMinMax,
minPSub,
TEnd,
...)
}
\arguments{
  \item{aTwDEMC}{former run of twDEMCBlockInt}
  \item{qPop}{numeric vector (nPop): probability of samples in subspace within entire space}
  \item{nGen}{the number of generations for twDEMCBlock}
  \item{controlTwDEMC}{list argument to \code{\link{twDEMCBlockInt}} containing entry thin}
  \item{m0}{minimum number of samples in step for extending runs}
  \item{nRowsMin}{number of rows in sub to sample at minimum, must contain enough samples for robust estimation of mean log-Density of subspace}
  \item{nRowsMinMax}{number of rows in sub to include as maximum mcNewMinN}
  \item{minPSub}{minimum proportion of a subsample in splitting´}
  \item{TEnd}{numeric vector (nResComp) target temperature}
  \item{\dots}{further arguments to \code{\link{twDEMCBlockInt}}, such as TEnd}
}
\details{A first estimate of the proportion of samples from different subspaces
are the initial percentiles qp. 
The proportion of the samples from different subspaces is estimated
by the proportions of integrated density of the subspaces.
These proportions are esimated by average density over all samples multiplied by an
estimate proportional to the volume: the initial quantiles.}
\value{For each population, a list with entries
\item{resTwDEMC}{the runs with new samples of class twDEMCPops}
\item{pSubs}{the quantiles of samples within given subspace}
\item{subPercChange}{numeric vector: ratio of estimated proportion in limiting distribution to proportion of initial proportion}
\item{resTwDEMCMinN}{twDEMCPops with last (nRowsMin <- max( 2*m0, nSampleMin )) rows}
\item{nPopsApp}{numeric vectgor (nSpace): number of spaces for which samples are appended}}

\author{Thomas Wutzler}




\examples{
data(den2dCorEx)
aTwDEMC <- den2dCorEx$mcSubspaces0
getSpacesPop(aTwDEMC)
getNSamples(aTwDEMC)

#mtrace(divideTwDEMCStep)
res <- divideTwDEMCStep(aTwDEMC, nGen=256, dInfos=list(list(fLogDen=den2dCor))
	,  debugSequential=TRUE
	#,  controlTwDEMC=list(DRGamma=0.05) 
)
getNSamples(res$resTwDEMC)

#windows(record=TRUE)
plot( as.mcmc.list(stackPopsInSpace(aTwDEMC)), smooth=FALSE)
#mc1 <- stackPopsInSpace(res$resTwDEMC, mergeMethod="stack")
mc1 <- stackPopsInSpace(res$resTwDEMC, mergeMethod="random")
plot( as.mcmc.list(mc1), smooth=FALSE) # note how the distribution shifted
#plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=4)), smooth=FALSE)
#plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=30)), smooth=FALSE)

ss <- stackChains(concatPops(mc1))			# the new sample
ss0 <- stackChains(concatPops(squeeze(stackPopsInSpace(aTwDEMC),length.out=getNSamples(mc1)[1]))) # initial sample of same length
plot( b ~ a, as.data.frame(ss0), xlim=c(-0.5,2), ylim=c(-20,40) )
plot( b ~ a, as.data.frame(ss), xlim=c(-0.5,2), ylim=c(-20,40) ) # not that more samples are in the region of interest

.tmp.f <- function(){
	#barplot(res[[1]]$subPercChange, names.arg=seq_along(res[[1]]$subPercChange) )
	#str(res[[1]]$lowerParBounds)
	ssImpPops1 <- ssImpPops <- abind( lapply( res, "[[", "sample"), rev.along=0 )
	plot(density( ssImpPops[,"a",1]));lines(density( ssImpPops[,"a",2]),col="green"); lines(density( ss0[,"a"]),col="blue")
	plot( b ~ a, as.data.frame(ss0), xlim=c(-0.5,2), ylim=c(-20,40) ); points(0.8, 0, col="red" )
	#plot( b ~ a, as.data.frame(ssImpPops[,,2]) ); points(xyMax[1], xyMax[2], col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,2]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(xyMax[1], xyMax[2], col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,1]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(xyMax[1], xyMax[2], col="red" )
	#plot(density( ssImpPops[,"b",1]));lines(density( ssImpPops[,"b",2]),col="green"); lines(density( ss0[,"b"]),col="blue")
	ssImpPops2 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops1[,,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE ), "[[", "sample"), rev.along=0 )
	#mtrace(divideTwDEMC)
	#ssImpPops2 <- ssImpPops <- divideTwDEMC(ssImpPops1[,-1,], nGen=100, fLogDen=den2dCor, attachDetails=TRUE )
	ssImpPops3 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops2[,,], nGen=500, fLogDen=den2dCor ), "[[", "sample"), rev.along=0 )
	ssImpPops4 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops3[,,], nGen=500, fLogDen=den2dCor ), "[[", "sample"), rev.along=0 )
	ssImpPops5 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops4[,,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE ), "[[", "sample"), rev.along=0 )
	ssImpPops6 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops5[,,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE ), "[[", "sample"), rev.along=0 )
	str(ssImpPops)
}
}
