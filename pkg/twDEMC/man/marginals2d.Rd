\name{marginals2d}
\alias{marginals2d}
\title{marginals2d}
\description{integrate over all others but 2 columns}
\usage{
marginals2d(x, vars = c(2, 3), intCol = 1, n = round(sqrt(nrow(x))/12), dimnames = .colNames(x, vars), fAgg = mean, ...)
}
\arguments{
  \item{x}{numeric matrix}
  \item{vars}{index (integer or string): column to define integration boxes}
  \item{intCol}{index (integer or string): column to average}
  \item{n}{numeric vector: number of groups per variable}
  \item{dimnames}{dimension names of the result array}
  \item{fAgg}{function applied over the subsets of x[,intCol]}
  \item{\dots}{further arguments to fAgg}
}

\value{numeric matrix with each item representing application of fAgg to x[,intCol]
and dimnames equal to means of the classes}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{marginals1d}}
\code{\link{twDEMCInt}}}
\examples{
res <- marginals2d(sample,vars=c("kY","kO"),n=10)
tmp.f <- function(){
	library(lattice)
	levelplot( value~kY*kO, data=melt(res), col.regions=rev(heat.colors(100)), xlab=names(dimnames(res))[1], ylab=names(dimnames(res))[2] )
}
}
