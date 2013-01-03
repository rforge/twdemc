\name{transOrigPopt.mcmc.list}
\alias{transOrigPopt.mcmc.list}
\title{transOrigPopt mcmc list}
\description{Applies \code{\link{transOrigPopt.default}} to each entry of \code{normopt}.}
\usage{\method{transOrigPopt}{mcmc.list}(normpopt, poptDistr = parDistr[colnames(normpopt[[1]]$parms), 
    "trans"], parDistr, ...)}
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
\examples{FALSE}
