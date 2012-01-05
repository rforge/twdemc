\name{divideTwDEMCSteps}
\alias{divideTwDEMCSteps}
\title{divideTwDEMCSteps}
\description{run twDEMCBlock on subspaces}
\usage{
divideTwDEMCSteps(aTwDEMC,
nGen = 512,
controlTwDEMC = list(),
doRecordProposals = FALSE,
...,
nGenBatch = 256,
dumpfileBasename = "recover",
critSpecVarRatio = 20,
minPSub = 0.1,
subPercChangeCrit = 1.6,
maxNSample = 128,
pThinPast = 0.5)
}
\arguments{
  \item{aTwDEMC}{former run of twDEMCBlockInt}
  \item{nGen}{the number of generations for twDEMCBlock}
  \item{controlTwDEMC}{list argument to \code{\link{twDEMCBlock}} containing entry thin}
  \item{doRecordProposals}{if TRUE then an array of each proposal together with the results of fLogDen are recorded and returned in component Y}
  \item{\dots}{further arguments to \code{\link{twDEMCBlock}}}
  \item{nGenBatch}{
}
  \item{dumpfileBasename}{
}
  \item{critSpecVarRatio}{
}
  \item{minPSub}{
}
  \item{subPercChangeCrit}{if all subPercChange of all sub-populations are below this value in two last batches, may assume convergence and skip further batches}
  \item{maxNSample}{if given a value, then looking for subspaces is done on a subsample of given length for efficiency (see \code{\link{getSubSpaces}})}
  \item{pThinPast}{in each step thin the past to given fraction before appending it}
}



\author{Thomas Wutzler}




\examples{
data(den2dCorEx)
#trace("twDEMCBlockInt", recover)
mc0 <- den2dCorEx$mcSubspaces0
#plot( as.mcmc.list(mc0), smooth=FALSE )
#plot( 
#mc0$pops[[ length(mc0$pops) ]]$T0 <- mc0$pops[[ length(mc0$pops) ]]$TEnd <- 100
mc0$blocks[[1]]$fUpdateBlock <- updateBlockTwDEMC	# if beeing traced
#mc0 <- res$resTwDEMC
#mtrace(divideTwDEMCSteps)
res <- divideTwDEMCSteps(mc0
	#, nGen=256*2+32
	, nGen=256*5
	, nGenBatch=256
	, dInfos=list(list(fLogDen=den2dCor))
	,  debugSequential=TRUE
	,  controlTwDEMC=list(DRgamma=0.1)
	, minPSub=0.05
)
getNSamples(res$resTwDEMC)
sort(res$subPercChange, decreasing=TRUE)
#str(controlTwDEMC)
#str(list(...)$controlTwDEMC)	
lapply(res$resTwDEMC$pops,"[[","splits")
#.getParBoundsPops(res$resTwDEMC$pops)

#windows(record=TRUE)
plot( as.mcmc.list(stackPopsInSpace(mc0)), smooth=FALSE)
#mc1 <- stackPopsInSpace(res$resTwDEMC, mergeMethod="stack")
mc1 <- stackPopsInSpace(res$resTwDEMC, mergeMethod="random")
plot( as.mcmc.list(mc1), smooth=FALSE) # note how the distribution shifted
#plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=4)), smooth=FALSE)
#plot( as.mcmc.list(concatPops(res$resTwDEMC,minPopLength=30)), smooth=FALSE)

ss <- stackChains(concatPops(mc1))			# the new sample
ss0 <- stackChains(concatPops(mc0))
.nr <- max( nrow(ss), nrow(ss0) )
ss <- ss[1:.nr,]
ss0 <- ss0[1:.nr,]
plot( b ~ a, as.data.frame(ss0), xlim=c(-0.5,2), ylim=c(-20,40) )
plot( b ~ a, as.data.frame(ss), xlim=c(-0.5,2), ylim=c(-20,40) ) # not that more samples are in the region of interest
points( 0.8,0.8, col="red" )
plot( b ~ a, as.data.frame(ss), xlim=c(-2,2), ylim=80*c(-20,40) ) # not that more samples are in the region of interest
points( 0.8,0.8, col="red" )

gridx <- a <- seq(-0.5,2,length.out=91)
gridy <- seq(-20,+40,length.out=91)
gridX <- expand.grid(gridx, gridy)
luden <- apply( gridX, 1, den2dCor ) 
mLuden <- matrix(luden,nrow=length(gridx))
image( gridx, gridy,  matrix(exp(luden),nrow=length(gridx)), col = rev(heat.colors(100)), xlab="a", ylab="b" )
points( b ~ a, as.data.frame(ss), xlim=c(-0.5,2), ylim=c(-20,40) ) # not that more samples are in the region of interest
points( 0.8,0.8, col="blue" )

}
