\name{ggplotChainPopMoves}
\alias{ggplotChainPopMoves}
\title{ggplotChainPopMoves}
\description{Plot boxplots of the distribution of first and last fifth of steps for each population }
\usage{
ggplotChainPopMoves(resB, iChains = rep(1:ncol(resB$temp), each = ncol(resB$rLogLik)\%/\%ncol(resB$temp)), doSort = TRUE)
}
\arguments{
  \item{resB}{the twDEMC to examine}
  \item{iChains}{mapping from chain to population}
  \item{doSort}{if TRUE result $rLogLik is sorted}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{plotMarginal2D}}
\code{\link{twDEMCInt}}}

