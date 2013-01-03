\name{transNormPopt.default}
\alias{transNormPopt.default}
\title{transNormPopt default}
\description{Transform vectors from original to normal scale.}
\usage{\method{transNormPopt}{default}(popt, poptDistr = parDistr[names(normpopt), 
    "trans"], parDistr, ...)}
\arguments{
  \item{popt}{numerical vector/array with values at untransformed scale}
  \item{poptDistr}{character vector/array of kind of transformation ("lognorm"/"logitnorm")
values with other characters indicate no transformation
positions must match the positions in normpopt}
  \item{parDistr}{alternative way of specifying poptDistr: 
dataframe with parameter names in column names and column trans, such as provided by \code{\link{twQuantiles2Coef}} }
  \item{\dots}{
}
}
\details{Argument \code{poptDistr} should have the same dimensions as normpopt. However, it is recycled
By this way it is possible to specify only one value, or vector corresponding to the rows of a matrix.}


\author{Thomas Wutzler}



\seealso{\code{\link{transOrigPopt.default}}
\code{\link{twDEMCBlockInt}}}
\examples{FALSE}
