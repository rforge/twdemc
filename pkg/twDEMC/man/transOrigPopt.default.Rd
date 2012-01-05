\name{transOrigPopt.default}
\alias{transOrigPopt.default}
\title{transOrigPopt default}
\description{Transform vectors from normal to original scale.}
\usage{
\method{transOrigPopt}{default}(normpopt,
poptDistr = eval(parse(text = "parDistr$trans[names(normpopt)]")),
...)
}
\arguments{
  \item{normpopt}{numerical vector/array with values at transformed, i.e. normal, scale}
  \item{poptDistr}{character vector/array of kind of transformation ("lognorm"/"logitnorm")
values with other characters indicate no transformation
default assumes vector \code{parDistr$trans} in environement }
  \item{\dots}{
}
}
\details{This generic method is provided for in several forms for first argument. \itemize{
\item{ as a named vector: this method  } 
\item{ as a matrix: \code{\link{transOrigPopt.matrix}}  } 
\item{ as a coda's mcmc.list: \code{\link{transOrigPopt.mcmc.list}}  }
\item{ as a list of class twDEMC \code{\link{transOrigPopt.twDEMC}}  }
}

There are further methods deal with parameter transformations. \itemize{
\item{ transforming from original to normal scale: \code{\link{transNormPopt.default}}  } 
\item{ calculating mu and sigma at normal scale from quantiles: \code{\link{twQuantiles2Coef}}  } 
\item{ constraing the result list of \code{\link{twQuantiles2Coef}}, and adding variance-covariance matrix: \code{\link{twConstrainPoptDistr}}  }
\item{ calculating the density based on the distribution arguments \code{\link{dDistr}} }
}

Argument \code{poptDistr} should have the same dimensions as normpopt. However, it is recycled.
By this way it is possible to specify only one value, or vector corresponding to the rows of a matrix.}
\value{Normpopt with some values transformed by exp (poptDist=="lognorm") or plogis (poptDistr=="logitnorm").}

\author{Thomas Wutzler}



\seealso{\code{\link{twVarDistrLevels}}
\code{\link{twDEMCBlockInt}}
\code{\link{transNormPopt.default}}}
\examples{
upperBoundProb = 0.99	# quantile of the upper boundary
parmsBounds = list(		# mode and upper bound
	A0 = c(10,15)		
	,D0 = c(10, 100)
	,C0 = c(0.6,0.8)
)
varDistr <- twVarDistrVec( names(parmsBounds) )	# by default assumed normal
varDistr["D0"] <- "lognorm"
varDistr["C0"] <- "logitnorm"
parDistr <- twQuantiles2Coef( parmsBounds, varDistr, upperBoundProb=upperBoundProb, useMedian=FALSE )
parDistr
all.equal( upperBoundProb, pnorm(parmsBounds$A0[2], parDistr$mu["A0"], parDistr$sigmaDiag["A0"] ) )

# transform entire parameter vectors between scales
pOrig <- transOrigPopt( parDistr$mu, parDistr$trans[names(parDistr$mu)] )
# note that transform of mu slighly differs from the mode for lognormal and logitnormal 
pOrig
# back-transform to normal scale
pBack <- transNormPopt( pOrig, parDistr$trans[names(pOrig)] )	
all.equal( parDistr$mu, pBack )

# plot quantiles for given distributions
pGrid <- seq(0.01,0.99,length.out=31)
plot( qnorm(pGrid, mean=parDistr$mu["D0"], sd=parDistr$sigmaDiag["D0"]) ~ pGrid)
plot( qlnorm(pGrid, mean=parDistr$mu["D0"], sd=parDistr$sigmaDiag["D0"]) ~ pGrid); abline(h=parmsBounds[["D0"]][1], col="grey")

# plot densities for D0 parameter ranges
dGrid <- seq(3, 80, length.out=100)
denOrig1 <- dlnorm(dGrid, mean=parDistr$mu["D0"], sd=parDistr$sigmaDiag["D0"]) 
plot( denOrig1 ~ dGrid, type="l"); abline(v=parmsBounds[["D0"]][1], col="grey")

# now plot the same using a grid on normal scale, transforming them to original scale
dNormGrid <- seq( parDistr$mu["D0"]-2*parDistr$sigmaDiag["D0"], parDistr$mu["D0"]+2*parDistr$sigmaDiag["D0"], length.out=100)
dOrigGrid <- transOrigPopt(dNormGrid, parDistr$trans["D0"])
all.equal( dNormGrid, transNormPopt(dOrigGrid, parDistr$trans["D0"]) )
denOrig2 <- dDistr( dOrigGrid, parDistr$mu["D0"],parDistr$sigmaDiag["D0"],trans=parDistr$trans["D0"] )
points( denOrig2 ~ dOrigGrid, col="blue" )
}
