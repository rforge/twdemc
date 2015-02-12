\name{initZtwDEMCSub.matrix}
\alias{initZtwDEMCSub.matrix}
\title{initZtwDEMCSub matrix}
\description{generates an appropriate initial sample of parameter vectors for twDEMC from subsampling a matrix}
\usage{\method{initZtwDEMCSub}{matrix}(x, ..., vars = colnames(x), nChainPop = 4, 
    nPop = 1, m0 = calcM0twDEMC(length(unique(vars)), nChainPop))}\arguments{
  \item{x}{the mcmc matrix to subsample (nCases x nParms)}
  \item{\dots}{
}
  \item{vars}{which variables to keep}
  \item{nChainPop}{number of chains per population}
  \item{nPop}{number of populations}
  \item{m0}{number of required cases for initialization}
}

\value{an array of dimension suitable for Zinit for twDEMCBlockInt}

\author{Thomas Wutzler}



\seealso{\code{\link{initZtwDEMCNormal}}}
\examples{
data(twLinreg1)
attach( twLinreg1 )

.nChainPop=4
.nPop=2
.nPar=length(theta0)
Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChainPop=.nChainPop, nPop=.nPop)
ss <- do.call( rbind, twMisc::twListArrDim( Zinit))	#stack the chains 
ZinitSub <- initZtwDEMCSub(ss, nChainPop=.nChainPop, nPop=.nPop)
head(Zinit[,,1])
all.equal(dim(Zinit),dim(ZinitSub))

detach()
}
