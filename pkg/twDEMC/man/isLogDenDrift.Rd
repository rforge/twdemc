\name{isLogDenDrift}
\alias{isLogDenDrift}
\title{isLogDenDrift}
\description{check whether first quartile all the logDensities is significantly smaller than last quartile }
\usage{isLogDenDrift(logDenT, dInfos, alpha = 0.05, maxDrift = 0.3)}
\arguments{
  \item{logDenT}{numeric array (nStep x nResComp): logDensity (highest are best)}
  \item{dInfos}{list of lists with entry resCompPos (integer vector) specifying the position of result components for each density}
  \item{alpha}{the significance level for a difference}
  \item{maxDrift}{difference in LogDensity, below which no drift is signalled}
}
\details{Because of large sample sizes, very small differences may be significantly different.
Use argument minDiff to specify below which difference a significant difference is not regarded as drift.}
\value{TRUE if any of the logDensities are a significantly greater in the fourth quantile compared to the first quantile of the samples}

\author{Thomas Wutzler}




\examples{
data(twdemcEx1)
logDenT <- calcTemperatedLogDen( twdemcEx1, getCurrentTemp(twdemcEx1))
isLogDenDrift(logDenT, twdemcEx1$dInfos )
}
