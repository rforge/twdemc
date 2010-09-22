\name{getDiffLogLik.twDEMCProps}
\alias{getDiffLogLik.twDEMCProps}
\title{getDiffLogLik twDEMCProps}
\description{Extract the Differences in Log-Likelihood between accepted states and proposals.}
\usage{
getDiffLogLik.twDEMCProps(Y, resCols, nLastSteps = 128, temp = 1, ...)
}
\arguments{
  \item{Y}{matrix of proposals with row "accepted" and first step (column) initial state, rows: results components of fLogLik and third dimension chains.}
  \item{resCols}{the rows of Y with result components, either names or positions}
  \item{nLastSteps}{number of last steps of Y for which to extract diffs}
  \item{temp}{numeric matrix (resComp x pop): the temperature applied to the difference of LogLiks}
  \item{\dots}{further arguments passed to \code{\link{getAcceptedPos.twDEMCProps}}}
}

\value{numeric array ( component x nLastSteps x chain ) of Lp-La, the L}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





