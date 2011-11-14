\name{stackChainsPop.twDEMC}
\alias{stackChainsPop.twDEMC}
\title{stackChainsPop twDEMC}
\description{Combine MarkovChains of each population of a twDEMC. }
\usage{
\method{stackChainsPop}{twDEMC}(x, ..., varInRows = FALSE)
}
\arguments{
  \item{x}{
}
  \item{\dots}{
}
  \item{varInRows}{set to TRUE if rows hold variables and columns steps (as in Zinit of twDEMC), defaults to variables in columns}
}

\value{Array with first column the logDensity rLogDen and the remaining columns the variables
, rows are steps, third dimension is the population 
(but see argument \code{varInRows}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{stackChains.twDEMC}}
\code{\link{subChains.twDEMC}}}
\examples{
data(twdemcEx1)
res <- stackChainsPop(twdemcEx1)
str(res)
(tmp1 <- head(res[,,1]))

res2 <- stackChainsPop(twdemcEx1, varInRows=TRUE)
str(res2)
(tmp2 <- head(t(res2[,,1])))

identical( tmp1, tmp2)
}
