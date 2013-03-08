\name{ggplotDensity.twDEMCPops}
\alias{ggplotDensity.twDEMCPops}
\title{ggplotDensity twDEMCPops}
\description{Plotting the densities for each parameter.}
\usage{ggplotDensity.twDEMCPops(res, poptDistr = NULL, pMin = 0.05, 
    doTransOrig = FALSE, doDispLogDen = TRUE, idInfo = names(res$dInfos)[1], 
    nSamplesPop = 500, popNames = as.character(seq_along(res$pops)), 
    popCols = scale_colour_hue(1:.nPop)$palette(.nPop), legendTitle = "Populations")}
\arguments{
  \item{res}{the twDEMCPops object whose densities to plot}
  \item{poptDistr}{parameter Distributions for the prior, usually \code{poptDistr <- \link{twConstrainPoptDistr}(poptNames,HamerParameterPriors$parDistr )}}
  \item{pMin}{if > 0, the results are constrained to quantiles of rLogDen>percMin. Can avoid extremes}
  \item{doTransOrig}{if TRUE, parameters are translated to original scale}
  \item{doDispLogDen}{include density of LogDensitys}
  \item{idInfo}{
}
  \item{nSamplesPop}{thin to about these number of samples within each population}
  \item{popNames}{character vector (nPop): names of the populations displayed in colour legend}
  \item{popCols}{colors of the populations}
  \item{legendTitle}{
}
}



\author{Thomas Wutzler}



\seealso{\code{\link{plotMarginal2D}}
\code{\link{twDEMCBlockInt}}}
\examples{
data(twdemcEx1)  # from package twDEMC
res <- twdemcEx1
ggplotDensity.twDEMCPops( res )
ggplotDensity.twDEMCPops( res, doDispLogDen=FALSE )
}
