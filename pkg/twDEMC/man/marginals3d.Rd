\name{marginals3d}
\alias{marginals3d}
\title{marginals3d}
\description{integrate over all others but 3 columns}
\usage{marginals3d(x, vars = c(2, 3, 4), intCol = 1, n = round(nrow(x)^(1/3)/12), 
    dimnames = .colNames(x, vars), fAgg = mean, ...)}
\arguments{
  \item{x}{numeric matrix}
  \item{vars}{index (integer or string): column to define integration boxes}
  \item{intCol}{index (integer or string): column to average}
  \item{n}{numeric vector: number of groups per variable}
  \item{dimnames}{dimension names of the result array}
  \item{fAgg}{function applied over the subsets of x[,intCol]}
  \item{\dots}{further arguments to fAgg}
}

\value{numeric array with each item representing application of fAgg to x[,intCol]
and dimnames equal to means of the classes}

\author{Thomas Wutzler}



\seealso{\code{\link{marginals1d}}
\code{\link{twDEMCBlockInt}}}
\examples{
.tmp.f <- function(){
	res <- marginals3d(sample,vars=c("kY","kO","tLagLeaf"),n=8)
	#library(Rcmdr)
	ds <- melt(res)
	scatter3d(ds$kY, ds$tLagLeaf, ds$kO
		, surface=FALSE
		,bg="white", axis.scales=TRUE, grid=TRUE, ellipsoid=FALSE
		, xlab="kY",ylab="tLagLeaf", zlab="kO"
		, point.col=rev(heat.colors(100))[round(rescale(ds$value,to=c(1,100)))]
		,sphere.size=1.5
	)
}

}
