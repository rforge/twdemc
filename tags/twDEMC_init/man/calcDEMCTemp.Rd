\name{calcDEMCTemp}
\alias{calcDEMCTemp}
\title{calcDEMCTemp}
\description{Calculates the temperature for an exponential decrease from \code{T0} to \code{Tend} after \code{nGen} steps. 	}
\usage{
calcDEMCTemp(T0, Tend = 1, nGen, iGen = 1:nGen)
}
\arguments{
  \item{T0}{the initial temperature (before the first step at iGen=0)}
  \item{Tend}{the temperature at the last step}
  \item{nGen}{the number of genrations	}
  \item{iGen}{the steps for which to calculate the Temperature	}
}

\value{vector of Temperatures corresponding to steps iGen}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twDEMCInt}}}
\examples{
plot( 1:100, calcDEMCTemp(T0=100,Tend=5,nGen=100) )	
}
