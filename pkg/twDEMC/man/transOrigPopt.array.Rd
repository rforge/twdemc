\name{transOrigPopt.array}
\alias{transOrigPopt.array}
\title{transOrigPopt array}
\description{Applies \code{\link{transOrigPopt.default}} to each column, i.e. variable, of \code{normopt}.}
\usage{\method{transOrigPopt}{array}(normpopt, poptDistr = parDistr[colnames(normpopt), 
    "trans"], parDistr, ...)}\arguments{
  \item{normpopt}{numerical array with values at transformed, i.e. normal, scale}
  \item{poptDistr}{character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt}
  \item{parDistr}{Alternative way of specifying poptDistr: 
dataframe with parameter names in column names and column trans, such as provided by \code{\link{twQuantiles2Coef}} }
  \item{\dots}{
}
}
\details{either poptDistr has names for each column name
or poptDistr has the same length as colnames(normpopt)}


\author{Thomas Wutzler}





