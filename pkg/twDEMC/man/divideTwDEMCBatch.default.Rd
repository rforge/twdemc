\name{divideTwDEMCBatch.default}
\alias{divideTwDEMCBatch.default}
\title{divideTwDEMCBatch default}
\description{iteratively run batches using divideTwDEMCBatch}
\usage{
\method{divideTwDEMCBatch}{default}(x, ..., nGen = 2 * 512, nGenBatch = 512, thinPastFac = 0.2, subPercChangeCrit = 1.6, argsFSplitPop = vector("list", dim(x)[3]))
}
\arguments{
  \item{x}{numeric array: rows:steps, col:parameters including logLik, 3: independent populations}
  \item{\dots}{further arguments to \code{\link{divideTwDEMC}}}
  \item{nGen}{number of generations}
  \item{nGenBatch}{number of generations within one batch}
  \item{thinPastFac}{thinning the past between batches to speed up localization between 0 (no past) and 1 (keep entire past)}
  \item{subPercChangeCrit}{if all subPercChange of all sub-populations are belwo this value in two last batches, may assume convergence and skip further batches}
  \item{argsFSplitPop}{for each population: list of arguments  passed \code{\link{getSubSpaces}} and further to \code{\link{findSplit}}, e.g. for passing order of variables to check in \code{iVars} and \code{jVarsVar}}
}

\value{A list with components 
\item{sample}{numeric array same as argument \code{x} but thinned and new sample appended}
\item{subPercChange}{matrix of numeric vectors giving the weights of sub-populations}
\item{maxSubPercChange}{numeric vector: calculated ratio of subspace weights, see argument \code{wSubFacMax}}
\item{resDivideTwDEMC}{result of last call to \code{\link{divideTwDEMC}}}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(den2dCorTwDEMC)
aTwDEMC <- 	thin(den2dCorTwDEMC, start=300)
aSample <- stackChainsPop(aTwDEMC)
#mtrace(divideTwDEMCBatch.default)
ssRes1 <- ssRes <- divideTwDEMCBatch(aSample, nGen=512*3, fLogDen=den2dCor )
#ssRes <- divideTwDEMCBatch(ssRes1, nGen=512*1, fLogDen=den2dCor )	#calling divideTwDEMCBatch.divideTwDEMCBatch
#mtrace(getSubSpaces)
#tmp <- getSubSpaces(ssRes1$sample[,-1,1])
#mtrace(divideTwDEMCBatch.divideTwDEMCBatch)
ssRes2 <- ssRes <- divideTwDEMCBatch(ssRes1, nGen=512*6, fLogDen=den2dCor )	#calling divideTwDEMCBatch.divideTwDEMCBatch

#tmp <- divideTwDEMC( ssRes1$sample[,,2,drop=FALSE], nGen=512, fLogDen=den2dCor)
#ssImpPops <- abind( tmp[[1]]$sample, tmp[[1]]$sample, rev.along=0) 

#wSubsB <- ssRes$wSubs[,1]
#wSubsB <- ssRes$wSubs[,2]
ssRes$maxSubPercChange
wSubsB <- ssRes$wSubs
sapply( lapply(wSubsB, quantile, probs = c(1, 0.25) ), function(entry){ entry[1] / entry[2] }) 

ssImpPops <- ssRes$sample
#ssImpPops <- ssRes1$sample
plot(density( ssImpPops[,"a",1]));lines(density( ssImpPops[,"a",2]),col="green"); lines(density( ss0[,"a"]),col="blue")
#ssImpPops <- ssRes1$sample
plot( b ~ a, as.data.frame(ssImpPops[,,1]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(0.8,0, col="red" )
plot( b ~ a, as.data.frame(ssImpPops[,,2]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(0.8,0, col="red" )
plot( b ~ a, as.data.frame(ssImpPops[,,1])); points(0.8,0, col="red" )
plot( b ~ a, as.data.frame(ssImpPops[,,2])); points(0.8,0, col="red" )
str(ssImpPops)
}
