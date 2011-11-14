\name{normConstLaplaceBridge}
\alias{normConstLaplaceBridge}
\title{normConstLaplaceBridge}
\description{Laplace bridge Approximation of the normalizing constant}
\usage{
normConstLaplaceBridge(sample, logPostSample, fLogPost, alpha = 0.05, nd = cov.rob(sample), M = 10000, logCest = normConstLaplace(sample, fLogPost, alpha, nd)["logCLV"])
}
\arguments{
  \item{sample}{MCMC sample of parameters, each row of the matrix is a parameter vector}
  \item{logPostSample}{the unnormalized posteror density of the sample}
  \item{fLogPost}{function to evaluate unnormalized Log-posterior for a given parameter vector}
  \item{alpha}{quantile defining the multivariate normal ellipsoid for volume correction}
  \item{nd}{normal density approximation: list with entries center and cov}
  \item{M}{number of additional simulations of posterior}
  \item{logCest}{initial estimate of the normalizing constant}
}

\value{log of the normalization constant}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
mu = c(0,0)
sd = c(1,2)
corr = diag(nrow=2)
corr[1,2] <- corr[2,1] <- 0.49
Sigma = diag(sd, nrow=length(sd)) \%*\% corr \%*\% diag(sd,nrow=length(sd)) 
n <- 1000
cTrue <- 10
sample <- rmvnorm(n,mean=mu,sigma=Sigma)
fLogPost <- function(theta){ log(cTrue) + dmvnorm(theta,mean=mu,sigma=Sigma, log=TRUE)}
logPostSample <- fLogPost(sample)
(CMW <- exp(normConstLaplaceBridge( sample, logPostSample, logPostSample )) )

#normd <- cov.rob(sample)
#(res <-  boot( sample, function(sample, i){normConstLaplace(sample[i,], fLogDen=fLogDen, nd=normd)}, 100))
#(res2 <-  boot( sample, function(sample, i){normConstLaplace(sample[i,], fLogDen=fLogDen)}, 100))

#normd <- cov.rob(sample)
#(res <-  boot( sample, function(sample, i){normConstLaplaceBridge(sample[i,], fLogDen=fLogDen, logDenSample, nd=normd)}, 100))
#(res2 <-  boot( sample, function(sample, i){normConstLaplaceBridge(sample[i,], fLogDen=fLogDen, logDenSample)}, 100))
}
