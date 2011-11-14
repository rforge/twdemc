\name{getDiffLogDen.twDEMCProps}
\alias{getDiffLogDen.twDEMCProps}
\title{getDiffLogDen twDEMCProps}
\description{Extract the Differences in LogDensity between accepted states and proposals.}
\usage{
getDiffLogDen.twDEMCProps(Y, resCols, nLastSteps = 128, temp = 1, ...)
}
\arguments{
  \item{Y}{matrix of proposals with row "accepted" and first step (column) initial state, rows: results components of fLogDen and third dimension chains.}
  \item{resCols}{the rows of Y with result components, either names or positions}
  \item{nLastSteps}{number of last steps of Y for which to extract diffs}
  \item{temp}{numeric matrix (resComp x pop): the temperature applied to the difference of LogDens}
  \item{\dots}{further arguments passed to \code{\link{getAcceptedPos.twDEMCProps}}}
}

\value{numeric array ( component x nLastSteps x chain ) of Lp-La, the L}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





