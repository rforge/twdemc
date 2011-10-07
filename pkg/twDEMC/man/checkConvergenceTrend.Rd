\name{checkConvergenceTrend}
\alias{checkConvergenceTrend}
\title{checkConvergenceTrend}
\description{checks whether the first and last fifth mean of populations differ significantly }
\usage{
checkConvergenceTrend(resB, iChains = rep(1:ncol(resB$temp), each = ncol(resB$rLogLik)\%/\%ncol(resB$temp)))
}
\arguments{
  \item{resB}{the twDEMC to examine}
  \item{iChains}{
}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(twdemcEx1)
# p.value for differenc in means for each population
(res <- checkConvergenceTrend(twdemcEx1))
# none of them has converged on a 5\% level
res < 0.05
#discard burnin 
(res <- checkConvergenceTrend(thin(twdemcEx1,start=40)))
# none of them has converged on a 5\% level
res < 0.05
}
