\name{twCoefLnormCi}
\alias{twCoefLnormCi}
\title{twCoefLnormCi}
\description{Calculates mu and sigma of the lognormal distribution from lower and upper quantile, i.e. conidence interval.}
\usage{twCoefLnormCi(lower, upper, sigmaFac = qnorm(0.99), isTransScale = FALSE)}
\arguments{
  \item{lower}{value at the lower quantile, i.e. practical minimum}
  \item{upper}{value at the upper quantile, i.e. practical maximum}
  \item{sigmaFac}{sigmaFac=2 is 95\% sigmaFac=2.6 is 99\% interval}
  \item{isTransScale}{if true lower and upper are already on log scale}
}

\value{named numeric vector: mu and sigma parameter of the lognormal distribution.}

\author{Thomas Wutzler}



\seealso{\code{\link{twQuantiles2Coef}}
\code{\link{transOrigPopt.default}}}
\examples{
mu=2
sd=c(1,0.8)
p=0.99
lower <- l <- qlnorm(1-p, mu, sd )		# p-confidence interval
upper <- u <- qlnorm(p, mu, sd )		# p-confidence interval
cf <- twCoefLnormCi(lower,upper)	
all.equal( cf[,"mu"] , c(mu,mu) )
all.equal( cf[,"sigma"] , sd )
}
