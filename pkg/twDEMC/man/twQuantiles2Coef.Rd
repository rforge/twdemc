\name{twQuantiles2Coef}
\alias{twQuantiles2Coef}
\title{twQuantiles2Coef}
\description{Calculating coefficients of transformed normal distributions from quantiles.}
\usage{
twQuantiles2Coef(parmsBounds,
varDistr,
upperBoundProb = 0.99,
useMedian = FALSE,
useCi = FALSE)
}
\arguments{
  \item{parmsBounds}{list parameters, each entry a numeric vector of length 2 specifying mode/median/lower quantile and upper quantile value}
  \item{varDistr}{character vector identifying the distribution, i.e. transformation to normal, for each parameter}
  \item{upperBoundProb}{probability for upper quantile in parmsBounds}
  \item{useMedian}{if TRUE, the first entry of parmsBounds specifies the median, insted of the mode}
  \item{useCi}{if TRUE, the first entry of parmsBounds specifies the lower quantile, insted of the mode}
}

\value{parameter distribution information, list entries \describe{
\item{trans}{character vector: type of distribtution (norm,lognorm,logitnorm)}
\item{mu}{numeric vector: distribution parameter mu, i.e. expected values at normal scale}
\item{sigmaDiag}{numeric vector: standard deviation for each parameter, i.e. sqrt(diagonal of matrix parameter sigma) multivariate distrubtion without correlations.}
\%\item{invsigma}{numeric matrix: inverse of the distribution parameter sigma on multivariate (transformed) normal distribution}
}}

\author{Thomas Wutzler}



\seealso{\code{\link{twCoefLnorm}}
\code{\link{twCoefLnormCi}}
\code{\link{twCoefLnormMLE}}
\code{\link{transOrigPopt.default}}
\code{\link{twDEMCBlockInt}}}

