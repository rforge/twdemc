\name{getCurrentBaseTemp.twDEMCPops}
\alias{getCurrentBaseTemp.twDEMCPops}
\title{getCurrentBaseTemp twDEMCPops}
\description{Get the Base Temperature, i.e. cost reduction factor at aggregated level}
\usage{\method{getCurrentBaseTemp}{twDEMCPops}(x,
 nObs = x$args$nObs, TFix = x$args$ctrlT$TFix, 
    ...)}\arguments{
  \item{x}{object of class twDEMCPops}
  \item{nObs}{number of observations for each result component, see \code{\link{calcBaseTemp}}}
  \item{TFix}{fixed temperatures for several result components, see \code{\link{calcBaseTemp}}}
  \item{\dots}{
}
}

\value{numeric scalar: Base Temperature at aggregated level}

\author{Thomas Wutzler}



\seealso{\code{\link{getCurrentTemp.twDEMCPops}}
\code{\link{calcBaseTemp}}
\code{\link{getNGen.twDEMCPops}}
\code{\link{subset.twDEMCPops}}
,\code{\link{twDEMCBlockInt}}}

