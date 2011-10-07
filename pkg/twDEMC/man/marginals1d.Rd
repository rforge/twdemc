\name{marginals1d}
\alias{marginals1d}
\title{marginals1d}
\description{integrate over all others but 1 column}
\usage{
marginals1d(x, vars = 2, intCol = 1, n = round(nrow(x)/12), dimnames = .colNames(x, vars), fAgg = mean, ...)
}
\arguments{
  \item{x}{numeric matrix}
  \item{vars}{index (integer or string): column to define integration boxes}
  \item{intCol}{index (integer or string): column to average}
  \item{n}{number of groups, defaults to about 10 observations per group}
  \item{dimnames}{dimension names of the result array}
  \item{fAgg}{function applied over the subsets of x[,intCol]}
  \item{\dots}{further arguments to fAgg}
}
\details{There are several methods to get calculate marginal densities for multivariate sample. \itemize{
\item{ integrate over all others but 1 column: this method  } 
\item{ integrate over all others but 2 column: \code{\link{marginals2d}}  } 
\item{ integrate over all others but 3 column: \code{\link{marginals3d}}  } 
\item{ plotting a 2D marginal distribution: \code{\link{plotMarginal2D}}  } 
\item{ plotting a 2D kernel density estimate: \code{\link{plot2DKDMarginals}}  } 
}}
\value{numeric vector with each item representing application of fAgg to x[,intCol]
and dimnames equal to means of the classes}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twDEMCInt}}}
\examples{
res <- marginals1d(sample,vars="kO",n=20)
plot( res ~ as.numeric(names(res)) )
}
