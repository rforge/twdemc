\name{ggplotDensity.twDEMC}
\alias{ggplotDensity.twDEMC}
\title{ggplotDensity twDEMC}
\description{Plotting the densities for each parameter.}
\usage{ggplotDensity.twDEMC(res, poptDistr = NULL, pMin = 0.05, doTransOrig = FALSE, 
    doDispLogDen = TRUE, idInfo = names(res$dInfos)[1])}
\arguments{
  \item{res}{the twDEMC whose densities to plot}
  \item{poptDistr}{parameter Distributions for the prior, usually \code{poptDistr <- \link{twConstrainPoptDistr}(poptNames,HamerParameterPriors$parDistr )}}
  \item{pMin}{if > 0, the results are constrained to quantiles of rLogDen>percMin. Can avoid extremes}
  \item{doTransOrig}{if TRUE, parameters are translated to original scale}
  \item{doDispLogDen}{include density of LogDensitys}
  \item{idInfo}{vector of column (names or indices) used to constrain the sample (logDen > pMin), defaults to first column}
}



\author{Thomas Wutzler}



\seealso{\code{\link{plotMarginal2D}}
\code{\link{twDEMCBlockInt}}}
\examples{
data(twdemcEx1)  # from package twDEMC
res <- concatPops(twdemcEx1)
ggplotDensity.twDEMC( res )
}
