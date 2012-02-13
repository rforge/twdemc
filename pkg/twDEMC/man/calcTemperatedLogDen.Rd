\name{calcTemperatedLogDen}
\alias{calcTemperatedLogDen}
\title{calcTemperatedLogDen}
\description{Generic method behind \code{\link{calcTemperatedLogDen.default}}}
\usage{calcTemperatedLogDen(...)}
\arguments{
  \item{\dots}{
}
}



\author{Thomas Wutzler}




\examples{
data(twdemcEx1)
logDen0 <- stackChains(concatPops(twdemcEx1)$resLogDen)
maxLogDen <- apply(logDen0,2,max)
temp = c(obs=20, parms=1)
logDen <- t( (t(logDen0)-maxLogDen)*temp +maxLogDen )

logDenT <- calcTemperatedLogDen(logDen,temp)
.exp <- logDen0
names(dimnames(.exp)) <- NULL 
names(dimnames(logDenT)) <- NULL 
all.equal( .exp, logDenT )
}
