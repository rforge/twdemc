\name{transOrigPopt.matrix}
\alias{transOrigPopt.matrix}
\title{transOrigPopt matrix}
\description{Applies \code{\link{transOrigPopt.default}} to each column of \code{normopt}.}
\usage{
\method{transOrigPopt}{matrix}(normpopt,
poptDistr = parDistr$trans[colnames(normpopt)],
parDistr,
...)
}
\arguments{
  \item{normpopt}{numerical matrx with values at transformed, i.e. normal, scale}
  \item{poptDistr}{character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt
Positions must match the positions in normpopt.
If given a single value, it is repeated.}
  \item{parDistr}{Alternative way of specifying poptDistr: 
dataframe with parameter names in column names and column trans, such as provided by \code{\link{twQuantiles2Coef}} }
  \item{\dots}{
}
}
\details{either poptDistr has names for each column name
or poptDistr has the same length as colnames(normpopt)}


\author{Thomas Wutzler}





