\name{devide.twDEMC}
\alias{devide.twDEMC}
\title{devide twDEMC}
\description{Devide a twDEMC result into to two subPopulations}
\usage{
\method{devide}{twDEMC}(x, split, ...)
}
\arguments{
  \item{x}{
}
  \item{split}{named numeric scalar: value and name of the parameter used to devide the sample}
  \item{\dots}{
}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{stackChains.twDEMC}}}
\examples{
data(twdemcEx1)
#x <- twdemcEx1
res <- devide(twdemcEx1,split=c(a=10.8))
str(res$right)
}
