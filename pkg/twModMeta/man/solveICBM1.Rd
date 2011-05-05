\name{solveICBM1}
\alias{solveICBM1}
\title{solveICBM1}
\description{solve the ODE of \code{\link{derivICBM1}}}
\usage{
solveICBM1(x0, times, parms, input, fFmAtmosphere = fmAtmosphere, modMeta = modMetaICBM1(), useRImpl = FALSE)
}
\arguments{
  \item{x0}{numeric vector or matrix at t=0}
  \item{times}{times at which explicit estimates for y are desired. The first value in times must be the initial time.}
  \item{parms}{list of model parameters}
  \item{input}{list with dataframes entries leaf and root each with columns yr and obs}
  \item{fFmAtmosphere}{
}
  \item{modMeta}{metaInformation from model. Pass for efficiency or when using different units.}
  \item{useRImpl}{flag indicating to use the R implementation instead of C implementation.}
}

\value{result of \code{\link{lsoda}}  }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link[deSolve]{lsoda}}}

