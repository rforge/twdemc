\name{dDistr}
\alias{dDistr}
\title{dDistr}
\description{Translate density from normal to original scale}
\usage{
dDistr(x, mu, sigma, trans)
}
\arguments{
  \item{x}{numeric vector of quantile at original scale for to calculate density}
  \item{mu}{numeric vector (recycled)}
  \item{sigma}{numeric vector (recycled)}
  \item{trans}{character vector: the Transformation to use (norm,lognorm,logitnorm)}
}
\details{To evaluate density at original, i.e. lognormal, logitnormal, scale
the density at transformed normal scale has to be multiplied with the 
Jacobian, i.e the derivative of the transformation}
\value{numeric vector of max length}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





