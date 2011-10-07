\name{mod2dCorExLog}
\alias{mod2dCorExLog}
\title{mod2dCorExLog}
\description{example 2d model function of normal density + sin}
\usage{
mod2dCorExLog(theta, xval, mu, Sigma, c = 500)
}
\arguments{
  \item{theta}{parameter vector with names x and y}
  \item{xval}{additional argument, passed by ... in logDenGaussian}
  \item{mu}{mean of the normal dist}
  \item{Sigma}{covariance matrix of the normal distribution}
  \item{c}{multiplicative constant}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twDEMCInt}}
\code{\link{logDenGaussian}}}
\examples{
mu = c(0,0)
sd = c(2,2)
corr = diag(nrow=2)
corr[1,2] <- corr[2,1] <- 0.9
Sigma = diag(sd, nrow=length(sd)) \%*\% corr \%*\% diag(sd,nrow=length(sd))


gridx <- gridy <- seq(-4,+4,length.out=61)
gridX <- expand.grid(gridx, gridy)
uden <- apply( gridX, 1, mod2dCorExLog, mu=mu, Sigma=Sigma ) 
image( gridx, gridy,  exp(matrix(uden,nrow=length(gridx))), col = rev(heat.colors(100)), xlab="x", ylab="y" )
image( gridx, gridy,  matrix(uden,nrow=length(gridx)), col = rev(heat.colors(100)), xlab="x", ylab="y" )

plot( density(uden) )
plot( sin(gridx) ~ gridx )
}
