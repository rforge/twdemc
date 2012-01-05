\name{normConstLaplace}
\alias{normConstLaplace}
\title{normConstLaplace}
\description{Laplace Approximation of the normalizing constant }
\usage{
normConstLaplace(sample,
fLogPost,
alpha = 0.05,
nd = cov.rob(sample))
}
\arguments{
  \item{sample}{MCMC sample of parameters, each row of the matrix is a parameter vector}
  \item{fLogPost}{function to evaluate log of unnormalized posterior density}
  \item{alpha}{quantile defining the multivariate normal ellipsoid for volume correction}
  \item{nd}{normal density approximation: list with entries center and cov}
}

\value{a list with estimates of the logs of the normalizing constant of the posterior
\item{logCL}{Laplace estimate}
\item{logCLV}{volume corrected Laplace estimate}}

\author{Thomas Wutzler}




\examples{
mu = c(0,0)
sd = c(1,2)
corr = diag(nrow=2)
corr[1,2] <- corr[2,1] <- 0.49
Sigma = diag(sd, nrow=length(sd)) \%*\% corr \%*\% diag(sd,nrow=length(sd)) 
n <- 1000
cTrue <- 10
sample <- rmvnorm(n,mean=mu,sigma=Sigma)
fPost <- function(theta){ cTrue* dmvnorm(theta,mean=mu,sigma=Sigma, log=FALSE)}
fLogPost <- function(theta){ log(cTrue) + dmvnorm(theta,mean=mu,sigma=Sigma, log=TRUE)}
exp(normConstLaplace( sample, fLogPost ))

#normd <- cov.rob(sample)
#(res <-  boot( sample, function(sample, i){normConstLaplace(sample[i,], fLogDen=fLogDen, nd=normd)}, 100))
#(res2 <-  boot( sample, function(sample, i){normConstLaplace(sample[i,], fLogDen=fLogDen)}, 100))
}
