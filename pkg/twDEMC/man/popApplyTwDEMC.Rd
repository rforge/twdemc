\name{popApplyTwDEMC}
\alias{popApplyTwDEMC}
\title{popApplyTwDEMC}
\description{Applying a function across all chains of one population for each case.}
\usage{
popApplyTwDEMC(x,
nPop,
FUN,
...)
}
\arguments{
  \item{x}{a matrix with columns chains or array with last dimension chain}
  \item{nPop}{number of populations}
  \item{FUN}{function to apply to population submatrix}
  \item{\dots}{further arguemtns to FUN}
}

\value{array with last dimenstion correponding to population}

\author{Thomas Wutzler}



\seealso{\code{\link{popMeansTwDEMC}}
\code{\link{subChains.twDEMC}}}
\examples{
data(twdemcEx1)
ex1c <- concatPops(twdemcEx1)
# mean logDen for each case, i.e. step, by population
nPop=getNPops(ex1c)
#mtrace(popApplyTwDEMC)
# applied to a matrix: pops in columns
popApplyTwDEMC( ex1c$logDen[,1,], nPop=nPop, apply, 1, mean )
# applied to a 3d array: pops in 3rd dimension
tmp <- popApplyTwDEMC( ex1c$logDen, nPop=nPop, apply, 1:2, mean )
str(tmp)

# stack logDen of block 1 for each population
tmp <- popApplyTwDEMC( ex1c$logDen[,1,], nPop=nPop, as.vector )
str(tmp)
#stack param columns by population
str(ex1c$parms)	#26 cases, 2 parameters, 2*4=8 chains
tmp <- popApplyTwDEMC( ex1c$parms, nPop=nPop, function(x){ abind(twListArrDim(x),along=1) })
all.equal( c(nrow(ex1c$parms)*4, 2,2), dim(tmp))
}
