\name{initZtwDEMCNormal}
\alias{initZtwDEMCNormal}
\title{initZtwDEMCNormal}
\description{Generate an initial population of states for \code{\link{twDEMCInt}}.}
\usage{
initZtwDEMCNormal(thetaPrior, covarTheta, nChains = 8, nPops = 2, m0 = calcM0twDEMC(length(thetaPrior), nPops, nChains), doIncludePrior = TRUE)
}
\arguments{
  \item{thetaPrior}{vector of parameters, point estimate}
  \item{covarTheta}{the a prior covariance of parameters }
  \item{nChains}{number of chains to run}
  \item{nPops}{number of independent populations among the chains }
  \item{m0}{number of initial states for each chain}
  \item{doIncludePrior}{If TRUE, then last sample of chain 1 to the prior estimate, 
which might be already a kind of best estimates by an optimization. }
}
\details{There are several methods to establish the initial population for a twDEMC run. \itemize{
\item{ drawing from a multivariate normal distribution: this method  } 
\item{ subsetting the result of a former twDEMC run: \code{\link{initZtwDEMCSub.twDEMC}}  } 
\item{ extending the result of a former twDEMC run to include more parameters: \code{\link{initZtwDEMCExt.twDEMC}}  } 
\item{ selecting the N closes points from a sequence of points in parameter space \code{\link{constrainNStack}}  } 
\item{ selecting the points inside a confindenc ellipsis in parameter \code{\link{constrainCfStack}}  } 
}}
\value{a matrix of number of parameters by number of individuals (d x m0 x Npop), with d dimension of theta}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twDEMCInt}}
\code{\link{calcM0twDEMC}}}

