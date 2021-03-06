\name{divideTwDEMC}
\alias{divideTwDEMC}
\title{divideTwDEMC}
\description{run twDEMC on subspaces}
\usage{
divideTwDEMC(aSample, nGen = 512, controlTwDEMC = list(), ..., minPSub = 0.05, isBreakEarly = FALSE, argsFSplitPop = vector("list", dim(aSample)[3]), nChainsPop = 4, m0 = calcM0twDEMC(ncol(aSample), nChains = nChainsPop), nrow0 = 4000, attachDetails = FALSE, nChainPar = 16, minNSamplesSub = 32, fCheckProblems = checkProblemsSpectral, argsFCheckProblems = list(), dumpfileBasename = "recover", subSpacesPop = list())
}
\arguments{
  \item{aSample}{numeric array: rows:steps, col:parameters including rLogLik in first column, 3: independent populations}
  \item{nGen}{the number of generations for twDEMC}
  \item{controlTwDEMC}{list argument to \code{\link{twDEMC}} containing entry thin}
  \item{\dots}{further arguments to \code{\link{twDEMC}}}
  \item{minPSub}{passed to \code{\link{getSubSpaces}}}
  \item{isBreakEarly}{passed to \code{\link{getSubSpaces}}}
  \item{argsFSplitPop}{for each population: list of arguments  passed \code{\link{getSubSpaces}} and further to \code{\link{findSplit}}, e.g. for passing order of variables to check in \code{iVars} and \code{jVarsVar}}
  \item{nChainsPop}{number of chains in subPopulations}
  \item{m0}{number of samples per chain to initialize subPopulations}
  \item{nrow0}{number of rows in initial overall sample}
  \item{attachDetails}{set TRUE to report upperParBounds, lowerParBounds, and pSubs per subPopulation}
  \item{nChainPar}{number of chains to run in parallel, good choice is 2*nCpu}
  \item{minNSamplesSub}{minimum number of records in a subspace sample, increase to avoid wrong estimation of weights}
  \item{fCheckProblems}{function applied to twDEMC results of each subspace}
  \item{argsFCheckProblems}{further arguments to argsFCheckProblems}
  \item{dumpfileBasename}{
}
  \item{subSpacesPop}{may provide previous results of \code{\link{getSubSpaces}} per population to save computing time}
}
\details{the first column of aSample records the logDensity of the sample for consitency. 
It is not used and may be initialized to any value. 

A first estimate of the proportion of samples from different subspaces
are the initial percentiles qp. 
The proportion of the samples from different subspaces is estimated
by the proportions of integrated density of the subspaces.
These proportions are esimated by average density over all samples multiplied by an
estimate proportional to the volume: the initial quantiles.}
\value{For each population, a list with entries
\item{sample}{numeric matrix: the sampled parameters (rows: cases, cols: 1: logDensity, others parameter dimensions}
\item{subPercChange}{numeric vector: estimated proportion in limiting distribution to proportion of initial proportion}
\item{upperParBounds}{list of named numeric vectors: the upper parameter bounds for the subspaces}
\item{lowerParBounds}{dito upperParBounds}
\item{iVars}{reordered index or parameter dimensions to check for split, see \code{\link{findSplit}}}
\item{jVarsVar}{reordered index or parameter dimensions to check for split, see \code{\link{findSplit}}}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(den2dCorTwDEMC)
aTwDEMC <- 	thin(den2dCorTwDEMC, start=300)
aSample <- stackChainsPop(aTwDEMC)
ss0 <- stackChains(aTwDEMC)
#aSamplePure <- aSample[,-1,]
#aSamplePurePop <- aSamplePure[,,1]
#res <- divideTwDEMC(aSample, nGen=100, fLogDen=den2dCor, attachDetails=TRUE )
#mtrace(divideTwDEMC)
res <- divideTwDEMC(aSample, nGen=500, fLogDen=den2dCor,  isBreakEarly=FALSE )
plot( b ~ a, as.data.frame(res[[1]]$sample), xlim=c(-0.5,2), ylim=c(-20,40) )
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
