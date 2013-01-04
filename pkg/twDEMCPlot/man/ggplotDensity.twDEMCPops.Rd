\name{ggplotDensity.twDEMCPops}
\alias{ggplotDensity.twDEMCPops}
\title{ggplotDensity twDEMCPops}
\description{Plotting the densities for each parameter.}
\usage{ggplotDensity.twDEMCPops(res, poptDistr = NULL, pMin = 0.05, 
    doTransOrig = FALSE, doDispLogDen = TRUE, idInfo = names(res$dInfos)[1], 
    nSamplesPop = 500)}
\arguments{
  \item{res}{the twDEMCPops object whose densities to plot}
  \item{poptDistr}{parameter Distributions for the prior, usually \code{poptDistr <- \link{twConstrainPoptDistr}(poptNames,HamerParameterPriors$parDistr )}}
  \item{pMin}{if > 0, the results are constrained to quantiles of rLogDen>percMin. Can avoid extremes}
  \item{doTransOrig}{if TRUE, parameters are translated to original scale}
  \item{doDispLogDen}{include density of LogDensitys}
  \item{idInfo}{
}
  \item{nSamplesPop}{thin to about these number of samples within each population}
}



\author{Thomas Wutzler}



\seealso{\code{\link{plotMarginal2D}}
\code{\link{twDEMCInt}}}
\examples{
data(twdemcEx1)  # from package twDEMC
res <- twdemcEx1
ggplotDensity.twDEMCPops( res )
ggplotDensity.twDEMCPops( res, doDispLogDen=FALSE )
}
