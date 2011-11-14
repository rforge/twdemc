\name{checkConvergenceTrend}
\alias{checkConvergenceTrend}
\title{checkConvergenceTrend}
\description{checks whether the first and last fifth mean of populations differ significantly }
\usage{
checkConvergenceTrend(resB, iChains = rep(1:ncol(resB$temp), each = ncol(resB$rLogDen)\%/\%ncol(resB$temp)))
}
\arguments{
  \item{resB}{the twDEMC to examine}
  \item{iChains}{
}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(twdemcEx1)
# p.value for difference in means for each population
(res <- checkConvergenceTrend(twdemcEx1))
# second has a trend on a 5\% level and did not yet converge
res > 0.05
}
