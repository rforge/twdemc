\name{dDistr}
\alias{dDistr}
\title{dDistr}
\description{Calculate density for transform of normal distribution.}
\usage{dDistr(x, mu, sigma, trans)}\arguments{
  \item{x}{numeric vector of quantile at original scale for to calculate density}
  \item{mu}{numeric vector (recycled)}
  \item{sigma}{numeric vector (recycled)}
  \item{trans}{factor  vector: the Transformation to use levels (norm,lognorm,logitnorm)}
}
\details{To evaluate density at original, i.e. lognormal, logitnormal, scale
the density at transformed normal scale has to be multiplied with the 
Jacobian, i.e the derivative of the transformation}
\value{numeric vector of length of maximum length of the arguments}

\author{Thomas Wutzler}



\seealso{\code{\link{transOrigPopt.default}}
\code{\link{twDEMCBlockInt}}}
\examples{
x <- seq(0.01,3,by=0.05)
mu=0
sigma=1
dx <- dDistr(x,mu,sigma,trans="lognorm")
plot( dx ~ x)
}
