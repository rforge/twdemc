\name{dinvchisq}
\alias{dinvchisq}
\title{dinvchisq}
\description{density function and for the (scaled) inverse-chi-squared distribution.}
\usage{dinvchisq(x, df, scale = 1/df, log = FALSE)}\arguments{
  \item{x}{vector of quantiles}
  \item{df}{degrees of freedom parameter, usually represented as nu}
  \item{scale}{scale parameter, usually represented as lambda.}
  \item{log}{Logical. If log=TRUE, then the logarithm of the density is returned.}
}



\author{Thomas Wutzler}



\seealso{\code{\link{dinvchisq}}}
\examples{
x <- dinvchisq(1,1,1)
x <- rinvchisq(10,1)

#Plot Probability Functions
x <- seq(from=0.1, to=5, by=0.01)
plot(x, dinvchisq(x,0.5,1), ylim=c(0,1), type="l", main="Probability Function",
        ylab="density", col="red")
lines(x, dinvchisq(x,1,1), type="l", col="green")
lines(x, dinvchisq(x,5,1), type="l", col="blue")
legend(3, 0.9, expression(paste(nu==0.5, ", ", lambda==1),
                paste(nu==1, ", ", lambda==1), paste(nu==5, ", ", lambda==1)),
        lty=c(1,1,1), col=c("red","green","blue"))	
}
