\name{initZtwDEMCNormal}
\alias{initZtwDEMCNormal}
\title{initZtwDEMCNormal}
\description{Generate an initial population of states for \code{\link{twDEMCBlockInt}}.}
\usage{initZtwDEMCNormal(thetaPrior, covarTheta, nChainPop = 4, nPop = 2, 
    m0 = ceiling(calcM0twDEMC(length(thetaPrior), nChainPop)/(m0FiniteFac)), 
    m0FiniteFac = 1, doIncludePrior = TRUE)}
\arguments{
  \item{thetaPrior}{numeric vector (nParm) of point estimate of the mean}
  \item{covarTheta}{numeric matrix (nParm x nParm) the covariance of parameters.
<< Alternatively, can be given as vector (nParm) of the diagonal, i.e. variances, if all dimensions are independent}
  \item{nChainPop}{number of chains to run}
  \item{nPop}{number of independent populations among the chains}
  \item{m0}{number of initial states for each chain}
  \item{m0FiniteFac}{use a factor smaller than 1 to increase default m0 to account for only a portion of proposal results in finite densities}
  \item{doIncludePrior}{If TRUE, then set last sample of chain 1 to the prior estimate, 
which might be already a kind of best estimates by an optimization. }
}
\details{There are several methods to establish the initial population for a twDEMC run. \itemize{
\item{ drawing from a multivariate normal distribution: this method  } 
\item{ subsetting the result of a former twDEMC run: \code{\link{initZtwDEMCSub.twDEMC}} or sample matrix \code{\link{initZtwDEMCSub.matrix}}  } 
\item{ extending the result of a former twDEMC run to include more parameters: \code{\link{initZtwDEMCExt.twDEMC}}  } 
\item{ selecting the N closes points from a sequence of points in parameter space \code{\link{constrainNStack}}  } 
\item{ selecting the points inside a confindenc ellipsis in parameter \code{\link{constrainCfStack}}  } 
\item{ replacing cases with in initial proposals that yield non-finite density \code{\link{replaceZinitNonFiniteLogDens}}  } 
\item{ general method for replacing cases in initial proposals \code{\link{replaceZinitCases}}  } 
}}
\value{a matrix of number of parameters by number of individuals (m0 x d x Npop), with d dimension of theta}

\author{Thomas Wutzler}



\seealso{\code{\link{twDEMCBlockInt}}
\code{\link{calcM0twDEMC}}}
\examples{
data(twLinreg1)
attach( twLinreg1 )

.nChainPop=4
.nPop=2
.nPar=length(theta0)
Zinit <- initZtwDEMCNormal( theta0, diag(sdTheta^2), nChainPop=.nChainPop, nPop=.nPop)
head(Zinit[,,1])
all.equal( c(calcM0twDEMC(.nPar,.nChainPop), .nPar, .nChainPop*.nPop), dim(Zinit) )

detach()
}
