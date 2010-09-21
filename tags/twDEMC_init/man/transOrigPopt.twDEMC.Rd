\name{transOrigPopt.twDEMC}
\alias{transOrigPopt.twDEMC}
\title{transOrigPopt twDEMC}
\description{Applies \code{\link{transOrigPopt.default}} to each column of parameters in \code{vtwdemc}.}
\usage{
transOrigPopt.twDEMC(vtwdemc, poptDistr = eval(parse(text = "parDistr$trans")), ...)
}
\arguments{
  \item{vtwdemc}{list of class twDEMC with $parms in transformed scale }
  \item{poptDistr}{character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt}
  \item{\dots}{
}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{transOrigPopt.default}}}

