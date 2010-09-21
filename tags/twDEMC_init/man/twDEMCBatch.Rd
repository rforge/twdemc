\name{twDEMCBatch}
\alias{twDEMCBatch}
\title{twDEMCBatch}
\description{Calls \code{\link{twDEMCBatchInt}} with arguments taken from attribute \code{batchCall} of \code{Zinit}.}
\usage{
twDEMCBatch(Zinit, ...)
}
\arguments{
  \item{Zinit}{the twDEMC object returned by \code{\link{twDEMCInt}}}
  \item{\dots}{Further arguments that are appended/overwrite entries of \code{batchCall}.
If Zinit is of class twDEMC, arguemnts  "resFLogLikX","logLikX" are set to NULL, so that these are taken from Zinit itself.}
}
\details{If twDEMC-list Zinit has no attribute batchCall, then \code{\link{twDEMCBatchInt}} is called with arguments provided.

If twDEMC-list Zinit has attribute \code{batchCall}, then it is re-executed with adjusted Zinit.
All further passed arguments will overwrite arguments in batchCall.
Hence it is possible to continue a run easily with \code{twDEMCBatch( resultPrev, nGen=nGenPrev+100 )}}
\value{result of \code{\link{twDEMCBatchInt}}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twDEMCBatchInt}}
\code{\link{twDEMCInt}}}

