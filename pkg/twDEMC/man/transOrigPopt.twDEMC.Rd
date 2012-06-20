\name{transOrigPopt.twDEMC}
\alias{transOrigPopt.twDEMC}
\title{transOrigPopt twDEMC}
\description{Applies \code{\link{transOrigPopt.default}} to each column of parameters in \code{vtwdemc}.}
\usage{
\method{transOrigPopt}{twDEMC}(normpopt,
poptDistr = parDistr$trans[colnames(normpopt$parms)],
parDistr,
...)
}
\arguments{
  \item{normpopt}{numerical matrx with values at transformed, i.e. normal, scale}
  \item{poptDistr}{character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt}
  \item{parDistr}{Alternative way of specifying poptDistr: 
dataframe with parameter names in column names and column trans, such as provided by \code{\link{twQuantiles2Coef}} }
  \item{\dots}{
}
}



\author{Thomas Wutzler}



\seealso{\code{\link{transOrigPopt.default}}}

