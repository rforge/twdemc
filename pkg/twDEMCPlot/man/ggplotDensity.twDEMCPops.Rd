\name{ggplotDensity.twDEMCPops}
\alias{ggplotDensity.twDEMCPops}
\title{ggplotDensity twDEMCPops}
\description{Plotting the densities for each parameter.}
\usage{ggplotDensity.twDEMCPops(res, poptDistr = NULL, pMin = 0.05, 
    doTransOrig = FALSE, doDispLogDen = TRUE, nSamplesPop = 500, 
    popNames = as.character(seq_along(res$pops)), popCols = scale_colour_hue(1:.nPop)$palette(.nPop), 
    legendTitle = "Populations", parNames = colnames(res$pops[[1]]$parms), 
    parmsBounds = NULL, parmsBoundsLt = "dashed", themeOpts = theme())}
\arguments{
  \item{res}{the twDEMCPops object whose densities to plot}
  \item{poptDistr}{parameter Distributions for the prior, usually \code{poptDistr <- \link{twConstrainPoptDistr}(poptNames,HamerParameterPriors$parDistr )}}
  \item{pMin}{if > 0, the results are constrained to quantiles of rLogDen>percMin. Can avoid extremes}
  \item{doTransOrig}{if TRUE, parameters are translated to original scale}
  \item{doDispLogDen}{include density of LogDensitys}
  \item{nSamplesPop}{thin to about these number of samples within each population}
  \item{popNames}{character vector (nPop): names of the populations displayed in colour legend}
  \item{popCols}{colors of the populations}
  \item{legendTitle}{
}
  \item{parNames}{names of the parameters to plot}
  \item{parmsBounds}{list <- colour -> numeric matrix (nParm x nLines) to be plottes as lines or vector (names or rownames must hold parameter names)}
  \item{parmsBoundsLt}{
}
  \item{themeOpts}{
}
}



\author{Thomas Wutzler}



\seealso{\code{\link{plotMarginal2D}}
\code{\link{twDEMCBlockInt}}}
\examples{
data(twdemcEx1)  # from package twDEMC
res <- twdemcEx1
ggplotDensity.twDEMCPops( res )
ggplotDensity.twDEMCPops( res, doDispLogDen=FALSE, popNames=c("scen 1","scen 2"), legendTitle="Scenarios" )
ggplotDensity.twDEMCPops( res, doDispLogDen=FALSE, parmsBounds=res$dInfos[[1]]$argsFLogDen$thetaPrior )  # indicate the prior mean
cols <- c("red","blue")
pBs <- structure( list(matrix(c(9,11,4.7, 5.2),byrow=TRUE,ncol=2,dimnames=list( c("a","b"),NULL))), names=cols[2])
ggplotDensity.twDEMCPops( res, doDispLogDen=FALSE, popCols=cols, parmsBounds=pBs )

#data(twTwoDenEx1); res <- subsetTail(int.ofMultiIntermediate())    # in runittwDEMC.R  # to test multiple cost 

}
