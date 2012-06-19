\name{initZtwDEMCExt.matrix}
\alias{initZtwDEMCExt.matrix}
\title{initZtwDEMCExt matrix}
\description{subsampling and extending twDEMC with new variables}
\usage{
\method{initZtwDEMCExt}{matrix}(Zinit1,
...,
thetaPrior,
covarTheta,
nChainPop = 4,
nPop = 1,
m0 = calcM0twDEMC(length(thetaPrior),
nChainPop))
}
\arguments{
  \item{Zinit1}{the matrix to subsample (nCases x nParms)}
  \item{\dots}{
}
  \item{thetaPrior}{numeric vector: mu of multivariate gaussian prior distribtuion}
  \item{covarTheta}{numeric vector: sigma of multivariate gaussian prior distrbituion}
  \item{nChainPop}{
}
  \item{nPop}{
}
  \item{m0}{number of required cases }
}
\details{the variables in thetaPrior that are part of vtwdemc are subsampled
the other variables are drawn from prior distribution
assuming no correlations between variables present and absent in vtwdemc and }


\author{Thomas Wutzler}



\seealso{\code{\link{initZtwDEMCNormal}}}

