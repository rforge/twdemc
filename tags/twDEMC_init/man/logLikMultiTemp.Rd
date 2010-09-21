\name{logLikMultiTemp}
\alias{logLikMultiTemp}
\title{logLikMultiTemp}
\description{LogLik function for the multi temperature test case }
\usage{
logLikMultiTemp(theta, logLikAccept = numeric(0), metropolisStepTemp = c(yPrior = 1), ..., thetaPrior = 0.8, sdTheta = 0.7, theta0 = 1, theta1 = 7, offset = -500, maxy = 12.964)
}
\arguments{
  \item{theta}{the parameter vector (here scalar).}
  \item{logLikAccept}{scalar: logLik for parms from revious run for two step Metropolis decision}
  \item{metropolisStepTemp}{numeric named vector: the temperature for internal metropolis step}
  \item{\dots}{any other arguments passed to fModel}
  \item{thetaPrior}{the prior estimate of the parameters}
  \item{sdTheta}{the standard deviation of the prior estimate}
  \item{theta0}{
}
  \item{theta1}{
}
  \item{offset}{
}
  \item{maxy}{
}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





