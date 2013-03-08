\name{ggplotDensity.poptDistr}
\alias{ggplotDensity.poptDistr}
\title{ggplotDensity poptDistr}
\description{Plotting the densities for each parameter.}
\usage{ggplotDensity.poptDistr(parDistr, pMin = 0.005, parmsBounds = NULL, 
    doTransOrig = TRUE)}
\arguments{
  \item{parDistr}{data.frame (nParm x c(mu, sigmaDiag)): parameter Distributions for the prior, usually \code{poptDistr <- twConstrainPoptDistr(poptNames,HamerParameterPriors$parDistr )}}
  \item{pMin}{range of the distribution from pMin to 1-pMin}
  \item{parmsBounds}{numeric vector (nParm x nLines) to be plottes as lines}
  \item{doTransOrig}{set to FALSE to display transform to normal scale}
}



\author{Thomas Wutzler}



\seealso{\code{\link{plotMarginal2D}}
\code{\link{twDEMCBlockInt}}}
\examples{
parmsBounds = do.call( rbind,list(		# mode and upper bound
	A0 = c(10,15)		
	,D0 = c(10, 100)
	,C0 = c(0.6,0.8)
))
colnames(parmsBounds) <- c("mode","upper")
varDistr <- twVarDistrVec( rownames(parmsBounds) )	# by default assumed normal
varDistr["D0"] <- "lognorm"
varDistr["C0"] <- "logitnorm"
parDistr <- twQuantiles2Coef( parmsBounds, varDistr, upperBoundProb=0.975 )

#mtrace(ggplotDensity.poptDistr)
ggplotDensity.poptDistr( parDistr, parmsBounds=parmsBounds )	
ggplotDensity.poptDistr( parDistr, parmsBounds=parmsBounds, doTransOrig=FALSE )	

}
