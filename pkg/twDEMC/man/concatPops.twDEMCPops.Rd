\name{concatPops.twDEMCPops}
\alias{concatPops.twDEMCPops}
\title{concatPops twDEMCPops}
\description{Concatenates all the chains of all subpopulations to one matrix of call \code{twDEMC}.}
\usage{
\method{concatPops}{twDEMCPops}(x,
...,
useThinning = TRUE,
minPopLength = NULL)
}
\arguments{
  \item{x}{
}
  \item{\dots}{
}
  \item{useThinning}{if TRUE thinning is used to make populations the same length, if FALSE they are cut to shortest population}
  \item{minPopLength}{integer scalar: if specified, populations with less samples than length.out are dropped}
}



\author{Thomas Wutzler}



\seealso{\code{\link{subset.twDEMCPops}}
,\code{\link{twDEMCBlockInt}}}

