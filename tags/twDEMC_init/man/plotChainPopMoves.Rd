\name{plotChainPopMoves}
\alias{plotChainPopMoves}
\title{plotChainPopMoves}
\description{Plot boxplots of the distribution of first and last fifth of steps for each population }
\usage{
plotChainPopMoves(resB, iChains = rep(1:ncol(resB$temp), each = ncol(resB$rLogLik)\%/\%ncol(resB$temp)), doSort = TRUE, ..., xlab = "Chains", ylab = "Log-Likelihood")
}
\arguments{
  \item{resB}{the twDEMC to examine}
  \item{iChains}{mapping from chain to population}
  \item{doSort}{if TRUE result $rLogLik is sorted}
  \item{\dots}{further arguements passed to matplot}
  \item{xlab}{
}
  \item{ylab}{
}
}
\details{\code{\link{ggplotChainPopMoves}} gives nicer results, but this functions i faster.}


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





