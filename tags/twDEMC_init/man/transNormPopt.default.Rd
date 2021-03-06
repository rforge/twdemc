\name{transNormPopt.default}
\alias{transNormPopt.default}
\title{transNormPopt default}
\description{Transform vectors from original to normal scale.}
\usage{
transNormPopt.default(popt, poptDistr = eval(parse(text = "parDistr$trans[names(normpopt)]")), ...)
}
\arguments{
  \item{popt}{numerical vector/array with values at untransformed scale}
  \item{poptDistr}{character vector/array of kind of transformation ("lognorm"/"logitnorm")
values with other characters indicate no transformation
default assumes vector parDistr$trans in environement }
  \item{\dots}{
}
}
\details{Argument \code{poptDistr} should have the same dimensions as normpopt. However, it is recycled
By this way it is possible to specify only one value, or vector corresponding to the rows of a matrix.}
\value{Argument \code{popt} with some values transformed by log (poptDist=="lognorm") or qlogis (poptDistr=="logitnorm").}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{transOrigPopt.default}}
\code{\link{twDEMCInt}}}

