\name{plotChainPopMoves}
\alias{plotChainPopMoves}
\title{plotChainPopMoves}
\description{Plot boxplots of the distribution of first and last fifth of steps for each population }
\usage{plotChainPopMoves(resB, iDen = 1, doSort = TRUE, ..., xlab = "Chains", 
    ylab = "LogDensity")}
\arguments{
  \item{resB}{the twDEMC to examine}
  \item{iDen}{index of the density component}
  \item{doSort}{if TRUE result sort $rLogDen in output}
  \item{\dots}{further arguements passed to matplot}
  \item{xlab}{
}
  \item{ylab}{
}
}
\details{\code{ggplotChainPopMoves} in twDEMCPlot gives nicer results, but this functions i faster.}


\author{Thomas Wutzler}



\seealso{\code{\link{plotThinned.mcmc.list}}
\code{\link{twDEMC}}}

